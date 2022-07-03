#pragma once

#include <algorithm>
#include <cassert>
#include <vector>

#include "../core/hexa_General.h"
#include "../math/hexa_Constants.h"
#include "hexa_Prewarpers.h"

namespace hexa
{
	enum class StateVariableType { HP, BP, BP1, LP, AP, LS, HS, tilt, BS };

	template <typename Type, typename Prewarper = TaylorPrewarper<Type>>
	class StateVariableFilter final
	{
		using FilterType = StateVariableType;
	public:
		StateVariableFilter() = default;

		//==============================================================================		
		void setCutoff(Type newCutoff) noexcept
		{
			if (utils::areSame(newCutoff, cutoff)) return;
			cutoff = std::clamp(newCutoff, Type(5.), Type(20.e3));
			update();
		}

		void setQ(Type newQ) noexcept
		{
			newQ = std::clamp(newQ, Type(0.001), Type(72));
			Type newR2 = 1 / newQ;
			if (utils::areSame(newR2, R2)) return;

			R2 = newR2;
			update();
		}

		void setBandWidth(Type newBW) noexcept
		{
			Type bw = std::clamp(newBW, Type(0.1), Type(3));

			Type bwM = std::exp2(bw / 2);
			Type newR2 = 2 * std::sinh(std::log2(pw.g(cutoff * bwM) / pw.g(cutoff / bwM)) * c<Type>::ln2 / 2);
			if (utils::areSame(newR2, R2)) return;

			R2 = newR2;
			update();
		}

		void setGain(Type newGain) noexcept
		{
			if (utils::areSame(newGain, gain)) return;
			gain = std::clamp(newGain, Type(-48.), Type(48.));
			update();
		}

		void setType(FilterType newType) noexcept
		{
			if (newType == type) return;
			type = newType;
			update();
		}

		//==============================================================================		
		Type getCutoff() const noexcept { return cutoff; }

		Type getR() const noexcept { return R2 / 2; }

		Type getGain() const noexcept { return gain; }

		Type getType() const noexcept { return type; }

		Type getSampleRate() const noexcept { return sampleRate; }


		//==============================================================================		
		void prepare(Type sRate, size_t numChannels, [[maybe_unused]] size_t maxBlockSize) noexcept
		{
			sampleRate = sRate;
			pw.setup(sampleRate);
	
			s1.resize(numChannels);
			s2.resize(numChannels);

			update();
			reset();
		}

		void process(const Type** inputs, Type** outputs, size_t nChans, size_t nFrames) noexcept
		{
			assert(nChans <= s1.size());
			assert(nChans <= s2.size());

			for (size_t ch = 0; ch < nChans; ++ch)
			{
				auto&& ls1 = s1[ch];
				auto&& ls2 = s2[ch];

				const Type* in = inputs[ch];
				Type* out = outputs[ch];

				for (size_t n = 0; n < nFrames; ++n)
				{
					out[n] = tick(in[n], ls1, ls2);
				}
			}
		}

		Type processSample(const Type& x, size_t ch)
		{
			assert(ch < s1.size());
			assert(ch < s2.size());
			return tick(x, s1[ch], s2[ch]);
		}

		// Same as processSample (introduced for brevity in complex processors)
		Type operator() (const Type& x, size_t ch)
		{
			assert(ch < s1.size());
			assert(ch < s2.size());
			return tick(x, s1[ch], s2[ch]);
		}

		void reset() noexcept
		{
			std::fill(s1.begin(), s1.end(), Type(0));
			std::fill(s2.begin(), s2.end(), Type(0));
		}

	private:
		Type tick(const Type& x, Type& s1, Type& s2)
		{
			Type b1 = g * x + s1, b2 = s2;

			// LU factorization and solution
			Type z2 = b2 - l21 * b1;
			Type u2 = z2 * u22Inv, u1 = (b1 - u12u22Inv * z2) * u11Inv;

			// Update integrators state
			s1 = 2 * u1 - s1;
			s2 = 2 * u2 - s2;

			return a1 * u1 + a2 * u2 + a0 * x;
		}

		//==============================================================================		
		void update() noexcept
		{
			Type m, m2, g1;
			switch (type)
			{
			case FilterType::LS:
				m = std::pow(Type(10), -gain / 80); m2 = m * m;
				g = pw.g(cutoff) * m;

				g1 = R2 * g + 1;
				l21 = -g / g1;
				u11Inv = 1 / g1;
				u22Inv = g1 / (g * (R2 + g) + 1);
				u12u22Inv = g * u22Inv;

				a1 = R2 / m2 - R2;
				a2 = 1 / (m2 * m2) - 1;
				a0 = 1;

				break;
			case FilterType::HS:
				m = std::pow(Type(10), gain / 80); m2 = m * m;
				g = pw.g(cutoff) * m;

				g1 = R2 * g + 1;
				l21 = -g / g1;
				u11Inv = 1 / g1;
				u22Inv = g1 / (g * (R2 + g) + 1);
				u12u22Inv = g * u22Inv;

				a1 = m2 * (1 - m2) * R2;
				a2 = 1 - m2 * m2;
				a0 = m2 * m2;

				break;
			case FilterType::tilt:
				m = std::pow(Type(10), gain / 40); m2 = m * m;
				g = pw.g(cutoff) * m;

				g1 = R2 * g + 1;
				l21 = -g / g1;
				u11Inv = 1 / g1;
				u22Inv = g1 / (g * (R2 + g) + 1);
				u12u22Inv = g * u22Inv;

				a1 = (1 - m2) * R2;
				a2 = 1 / m2 - m2;
				a0 = m2;

				break;
			case FilterType::BS:
				m = std::pow(Type(10), -gain / 40);
				g = pw.g(cutoff);

				g1 = R2 * g + 1;
				l21 = -g / g1;
				u11Inv = 1 / g1;
				u22Inv = g1 / (g * (R2 + g) + 1);
				u12u22Inv = g * u22Inv;

				a1 = (1 / m - m) * R2;
				a2 = 0;
				a0 = 1;

				break;
			case FilterType::HP:
				g = pw.g(cutoff);

				g1 = R2 * g + 1;
				l21 = -g / g1;
				u11Inv = 1 / g1;
				u22Inv = g1 / (g * (R2 + g) + 1);
				u12u22Inv = g * u22Inv;

				a1 = -R2;
				a2 = -1;
				a0 = 1;

				break;
			case FilterType::BP:
				g = pw.g(cutoff);

				g1 = R2 * g + 1;
				l21 = -g / g1;
				u11Inv = 1 / g1;
				u22Inv = g1 / (g * (R2 + g) + 1);
				u12u22Inv = g * u22Inv;

				a1 = 1;
				a2 = 0;
				a0 = 0;

				break;
			case FilterType::BP1:
				g = pw.g(cutoff);

				g1 = R2 * g + 1;
				l21 = -g / g1;
				u11Inv = 1 / g1;
				u22Inv = g1 / (g * (R2 + g) + 1);
				u12u22Inv = g * u22Inv;

				a1 = R2;
				a2 = 0;
				a0 = 0;

				break;
			case FilterType::LP:
				g = pw.g(cutoff);

				g1 = R2 * g + 1;
				l21 = -g / g1;
				u11Inv = 1 / g1;
				u22Inv = g1 / (g * (R2 + g) + 1);
				u12u22Inv = g * u22Inv;

				a1 = 0;
				a2 = 1;
				a0 = 0;

				break;
			case FilterType::AP:
				g = pw.g(cutoff);

				g1 = R2 * g + 1;
				l21 = -g / g1;
				u11Inv = 1 / g1;
				u22Inv = g1 / (g * (R2 + g) + 1);
				u12u22Inv = g * u22Inv;

				a1 = -2 * R2;
				a2 = 0;
				a0 = 1;

				break;
			}
		}

		//==============================================================================		
		Type sampleRate{ 44100. }, cutoff{ 440. }, gain{ 6. }, R2{ c<Type>::sqrt2 };
		FilterType type{ FilterType::LP };

		Type g{}, a1{}, a2{}, a0{};
		Type l21{}, u11Inv{}, u22Inv{}, u12u22Inv{};

		std::vector<Type> s1{ 2 }, s2{ 2 };

		Prewarper pw{};
	};
}
