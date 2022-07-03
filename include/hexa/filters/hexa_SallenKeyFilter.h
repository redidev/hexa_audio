#pragma once

#include <algorithm>
#include <cassert>
#include <vector>

#include "hexa_Prewarpers.h"

namespace hexa
{
	enum class SallenKeyFilterType { LP, HP, BP, BP1 };

	/**
	 * MIMO realization of a linear Sallen-Key filter
	 */
	template <typename Type, typename Prewarper = TaylorPrewarper<Type>>
	class SallenKeyFilter
	{
		using FilterType = SallenKeyFilterType;
	public:
		SallenKeyFilter() = default;

		//==============================================================================
		void setFrequency(Type frequencyHz) noexcept
		{
			cutoff = std::clamp(frequencyHz, Type(5), Type(20.e3));
			update<true, false, false>();
		}

		void setResonance(Type resonance) noexcept
		{
			reso = std::clamp(resonance, Type(0), Type(1));
			update<false, true, false>();
		}

		void setType(FilterType newType) noexcept
		{
			type = newType;
			update<true, true, true>();
		}

		//==============================================================================
		Type getCutoff() const noexcept { return cutoff; }

		Type getResonance() const noexcept { return reso; }

		Type getType() const noexcept { return type; }

		Type getSampleRate() const noexcept { return sampleRate; }

		//==============================================================================
		void prepare(Type sRate, size_t numChannels, [[maybe_unused]] size_t maxBlockSize) noexcept
		{
			sampleRate = sRate;
			pw.setup(sampleRate);

			st1.resize(numChannels);
			st2.resize(numChannels);

			update<true, true, true>();
			reset();
		}

		void process(const Type** inputs, Type** outputs, size_t nChans, size_t nFrames) noexcept
		{
			assert(nChans <= st1.size());
			assert(nChans <= st2.size());

			for (size_t ch = 0; ch < nChans; ++ch)
			{
				auto&& ls1 = st1[ch];
				auto&& ls2 = st2[ch];

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
			assert(ch < st1.size());
			assert(ch < st2.size());
			return tick(x, st1[ch], st2[ch]);
		}

		// Same as processSample (introduced for brevity in complex processors)
		Type operator() (const Type& x, size_t ch)
		{
			assert(ch < st1.size());
			assert(ch < st2.size());
			return tick(x, st1[ch], st2[ch]);
		}

		void reset() noexcept
		{
			std::fill(st1.begin(), st1.end(), Type(0));
			std::fill(st2.begin(), st2.end(), Type(0));
		}

	private:
		Type tick(const Type& x, Type& s1, Type& s2) noexcept
		{
			// RHS
			Type b1 = g * (a1 - a2) * x + s1;
			Type b2 = g * a2 * x + s2;

			// Calculate current state-space vector.
			Type u1 = m11 * b1 + m12 * b2;
			Type u2 = m21 * b1 + m22 * b2;

			// Update states.
			s1 = 2 * u1 - s1;
			s2 = 2 * u2 - s2;

			return c1 * u1 + c2 * u2 + c0 * x;
		}

		//==============================================================================
		template <bool updateFreq, bool updateReso, bool updateType>
		void update() noexcept
		{
			if constexpr (updateFreq) g = pw.g(cutoff);
			if constexpr (updateReso) k = 2 * reso;
			if constexpr (updateType)
			{
				a1 = type == FilterType::LP ? 1 : 0;
				a2 = type == FilterType::LP ? 0 : 1;
			}

			Type g1 = 1 + g;
			Type gk = g * k;
			Type d = 1 / (g1 * g1 - gk);

			// Calculate state-space matrix.
			m11 = d * (g1 - gk);
			m12 = d * -gk;
			m21 = d * g;
			m22 = d * g1;

			switch (type)
			{
			case FilterType::LP:	
			case FilterType::BP:
				c1 = 1;
				c2 = c0 = 0;

				break;
			case FilterType::HP:	
				c1 = 1;
				c2 = k - 1;
				c0 = a2;

				break;
			case FilterType::BP1: 
				c1 = 2 - k;
				c2 = c0 = 0;

				break;
			}
		}

		Type sampleRate{ 44100. }, cutoff{ 200. }, reso{ Type(0.95) };
		FilterType type{ FilterType::LP };

		Type g{}, k{}, a1{}, a2{};
		Type m11{}, m12{}, m21{}, m22{};
		Type c0{}, c1{}, c2{};
		std::vector<Type> st1{ 2 }, st2{ 2 };

		Prewarper<Type> pw;
	};
}
