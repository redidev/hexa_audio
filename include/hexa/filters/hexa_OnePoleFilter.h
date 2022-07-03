#pragma once

#include <algorithm>
#include <cassert>
#include <vector>

#include "../core/hexa_General.h"
#include "hexa_Prewarpers.h"

namespace hexa
{
	enum class OnePoleType { LP, HP, AP, LS, HS, tilt };

	template <typename Type, typename Prewarper = TaylorPrewarper<Type>>
	class OnePoleFilter final
	{
		using FilterType = OnePoleType;
	public:
		OnePoleFilter() = default;

		//==============================================================================		
		void setCutoff(Type newCutoff) noexcept
		{
			if (utils::areSame(cutoff, newCutoff)) return;
			cutoff =  std::clamp(newCutoff, Type(5.), Type(20.e3));
			update();
		}

		void setGain(Type newGain) noexcept
		{
			if (utils::areSame(gain, newGain)) return;
			gain = std::clamp(newGain, Type(-48.), Type(48.));
			update();
		}

		void setType(FilterType newType) noexcept
		{
			if (type == newType) return;
			type = newType;
			update();
		}

		//==============================================================================		
		Type getCutoff() const noexcept { return cutoff; }

		Type getGain() const noexcept { return gain; }

		Type getType() const noexcept { return type; }

		//==============================================================================		
		void prepare(Type sRate, size_t numChannels, [[maybe_unused]] size_t maxBlockSize) noexcept
		{
			sampleRate = sRate;
			pw.setup(sampleRate);
			s.resize(numChannels, Type());

			update();
			reset();
		}

		void process(const Type** inputs, Type** outputs, size_t nChans, size_t nFrames) noexcept
		{
			assert(nChans <= s.size());
	
			for (size_t ch = 0; ch < nChans; ++ch)
			{
				auto&& ls = s[ch];

				const Type* in = inputs[ch];
				Type* out = outputs[ch];

				for (size_t n = 0; n < nFrames; ++n)
				{
					out[n] = tick(in[n], ls);
				}
			}
		}

		Type processSample(const Type& x, size_t ch)
		{
			assert(ch < s.size());
			return tick(x, s[ch]);
		}

		// Same as processSample (introduced for brevity in complex processors)
		Type operator() (const Type& x, size_t ch)
		{
			assert(ch < s.size());
			return tick(x, s[ch]);
		}

		void reset() noexcept
		{
			std::fill(s.begin(), s.end(), Type(0));
		}

	private:
		Type tick(const Type& x, Type& ls)
		{
			auto v = G * (x - ls);
			auto u = v + ls;
			ls = u + v;

			return a1 * u + a0 * x;
		}

		//==============================================================================		
		void update() noexcept
		{
			Type m, m2;
			switch (type)
			{
			case FilterType::LP:
				g = pw.g(cutoff);
				G = g / (1 + g);

				a1 = 1;
				a0 = 0;

				break;
			case FilterType::HP:
				g = pw.g(cutoff);
				G = g / (1 + g);

				a1 = -1;
				a0 = 1;

				break;
			case FilterType::AP:
				g = pw.g(cutoff);
				G = g / (1 + g);

				a1 = 2;
				a0 = -1;

				break;
			case FilterType::HS:
				m = std::pow(Type(10), gain / 40); m2 = m * m;
				g = pw.g(cutoff) * m;
				G = g / (1 + g);

				a1 = 1 - m2;
				a0 = m2;

				break;
			case FilterType::tilt:
				m = std::pow(Type(10), gain / 20);
				g = pw.g(cutoff) * m;
				G = g / (1 + g);

				a1 = 1 / m - m;
				a0 = m;

				break;
			case FilterType::LS:
				m = std::pow(Type(10), -gain / 40); m2 = m * m;
				g = pw.g(cutoff) * m;
				G = g / (1 + g);

				a1 = 1 / m2 - 1;
				a0 = 1;

				break;
			}
		}

		//==============================================================================		
		Type sampleRate{ 44100. }, cutoff{ 440. }, gain{ +6. };
		FilterType type{ FilterType::LP };

		Type g{}, G{}, a1{}, a0{};
		std::vector<Type> s{};

		Prewarper pw{};
	};
}
