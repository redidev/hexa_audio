#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

#include "../core/hexa_General.h"

namespace hexa
{
	enum class RBJFilterType { LP, HP, BP, BP1, LS, HS, peak, notch, AP };

	/**
	 * Implementation of a classical bi-quad filter, based on the famous RBJ Cookbook paper
	 * by Robert Bristow-Johnson
	 */
	template <typename Type>
	class RBJFilter
	{
		using FilterType = RBJFilterType;

	public:
		RBJFilter() = default;

		//==============================================================================
		/** Sets the cutoff. */
		void setCutoff(Type freq)
		{
			if (utils::areSame(freq, cutoff)) return;
			cutoff = freq;
			update<true, false>();
		}

		/** Sets the filter quality. */
		void setQ(Type newQ)
		{
			Type newR = 1 / (newQ + newQ);
			if (utils::areSame(newR, R)) return;
			R = newR;
			update<false, false>();
		}

		void setGain(double gainDb)
		{
			if (utils::areSame(gainInDb, gainDb)) return;
			gainInDb = gainDb;
			update<true, true>();
		}

		/** Sets the type of the filter. */
		void setType(FilterType newType)
		{
			type = newType;
			update<true, true>();
		}

		//==============================================================================
		void prepare(Type sRate, size_t numChannels, [[maybe_unused]] size_t maxBlockSize) noexcept
		{
			sampleRate = sRate;

			st1.resize(numChannels);
			st2.resize(numChannels);

			update<true, true>();
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
		template <bool updateFreqParams, bool updateGainParams>
		void update() noexcept
		{
			if constexpr (updateGainParams)
			{
				ASqRt = std::pow(Type(10), gainInDb / 80);
				A = ASqRt * ASqRt;
			}

			if constexpr (updateFreqParams)
			{
				const auto w0 = c<Type>::twoPi * cutoff / sampleRate;
				cosw0 = std::cos(w0);
				sinw0 = std::sin(w0);
			}

			alpha = sinw0 * R;

			Type b0, b1, b2, a0, a1, a2;
			switch (type)
			{
			case FilterType::LP:
				b0 = (1 - cosw0) / 2;
				b1 = 1 - cosw0;
				b2 = (1 - cosw0) / 2;
				a0 = 1 + alpha;
				a1 = -2 * cosw0;
				a2 = 1 - alpha;
				break;
			case FilterType::BP:
				b0 = sinw0 / 2;
				b1 = 0;
				b2 = -sinw0 / 2.;
				a0 = 1 + alpha;
				a1 = -2 * cosw0;
				a2 = 1 - alpha;
				break;
			case FilterType::HP:
				b0 = (1 + cosw0) / 2;
				b1 = -(1 + cosw0);
				b2 = (1 + cosw0) / 2;
				a0 = 1 + alpha;
				a1 = -2 * cosw0;
				a2 = 1 - alpha;
				break;
			case FilterType::BP1:
				b0 = alpha;
				b1 = 0;
				b2 = -alpha;
				a0 = 1 + alpha;
				a1 = -2 * cosw0;
				a2 = 1 - alpha;
				break;
			case FilterType::LS:
				b0 = A * ((A + 1) - (A - 1) * cosw0 + 2 * ASqRt * alpha);
				b1 = 2 * A * ((A - 1) - (A + 1) * cosw0);
				b2 = A * ((A + 1) - (A - 1) * cosw0 - 2 * ASqRt * alpha);
				a0 = (A + 1) + (A - 1) * cosw0 + 2 * ASqRt * alpha;
				a1 = -2 * ((A - 1) + (A + 1) * cosw0);
				a2 = (A + 1) + (A - 1) * cosw0 - 2 * ASqRt * alpha;
				break;
			case FilterType::HS:
				b0 = A * ((A + 1) + (A - 1) * cosw0 + 2 * ASqRt * alpha);
				b1 = -2 * A * ((A - 1) + (A + 1) * cosw0);
				b2 = A * ((A + 1) + (A - 1) * cosw0 - 2 * ASqRt * alpha);
				a0 = (A + 1) - (A - 1) * cosw0 + 2 * ASqRt * alpha;
				a1 = 2 * ((A - 1) - (A + 1) * cosw0);
				a2 = (A + 1) - (A - 1) * cosw0 - 2 * ASqRt * alpha;
				break;
			case FilterType::peak:
				b0 = 1 + 1 * alpha * A;
				b1 = -2 * cosw0;
				b2 = 1 - 1 * alpha * A;
				a0 = 1 + alpha / A;
				a1 = -2 * cosw0;
				a2 = 1 - alpha / A;
				break;
			case FilterType::notch:
				b0 = 1;
				b1 = -2 * cosw0;
				b2 = 1;
				a0 = 1 + alpha;
				a1 = -2 * cosw0;
				a2 = 1 - alpha;
				break;
			case FilterType::AP:
				b0 = 1 - alpha;
				b1 = -2 * cosw0;
				b2 = 1 + alpha;
				a0 = 1 + alpha;
				a1 = -2 * cosw0;
				a2 = 1 - alpha;
				break;
			}

			// Calculate bi-quad coefficients
			a1Da0 = a1 / a0;
			a2Da0 = a2 / a0;
			b0Da0 = b0 / a0;
			b1Da0 = b1 / a0;
			b2Da0 = b2 / a0;
		}

		Type tick(const Type& x, Type& s1, Type& s2)
		{
			// Transposed canonical form (TDF-II).
			Type y = b0Da0 * x + s1;

			s1 = b1Da0 * x - a1Da0 * y + s2;
			s2 = b2Da0 * x - a2Da0 * y;

			return  y;
		}

		Type cutoff{ 200. }, gainInDb{ 6. }, R{ c<Type>::reciprSqrt2 }, sampleRate{ 44100. };
		FilterType type{ FilterType::LP };

		Type alpha, sinw0, cosw0, A, ASqRt;
		Type a1Da0{}, a2Da0{}, b0Da0{}, b1Da0{}, b2Da0{};

		//==============================================================================
		std::vector<Type> st1{ 2 }, st2{ 2 };
	};
}