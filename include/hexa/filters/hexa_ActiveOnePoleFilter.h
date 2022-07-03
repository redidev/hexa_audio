#pragma once

#include <vector>
#include <cmath>

#include "../math/hexa_Constants.h"

namespace hexa
{
	/**
	 * Active one pole filter with OTA (My challenge to Urs' one pole monster ;-) )
	 */
	template <typename Type>
	class ActiveOnePoleFilter
	{
	public:
		ActiveOnePoleFilter() = default;

		//==============================================================================
		void setFrequency(Type freqHz) noexcept
		{
			cutoff = std::clamp(freqHz, Type(20), Type(16.e3));
			update();
		}

		void setDrive(Type gainDb) noexcept
		{
			gain = std::pow(Type(10), gainDb / 20);
		}

		//==============================================================================
		Type getCutoff() const noexcept { return cutoff; }

		Type getDrive() const noexcept { return 20 * std::log10(gain); }

		Type getSampleRate() const noexcept { return sampleRate; }

		//==============================================================================
		void prepare(Type sRate, size_t numChannels, [[maybe_unused]] size_t maxBlockSize) noexcept
		{
			sampleRate = sRate;

			st.resize(numChannels);
			update();
			reset();
		}

		void process(const Type** inputs, Type** outputs, size_t nChans, size_t nFrames) noexcept
		{
			assert(nChans <= st.size());;

			for (size_t ch = 0; ch < nChans; ++ch)
			{
				auto&& ls = st[ch];

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
			assert(ch < st.size());
			return tick(x, st[ch]);
		}

		// Same as processSample (introduced for brevity in complex processors)
		Type operator() (const Type& x, size_t ch)
		{
			assert(ch < st.size());
			return tick(x, st[ch]);
		}

		void reset() noexcept
		{
			std::fill(st.begin(), st.end(), Type(0));
		}

	private:
		//==============================================================================
		Type tick(const Type& in, Type& s) noexcept
		{
			Type x = in * gain;

			// Set initial guess.
			Type yV = s;

			// Damped Newton (see Kelley monography)
			Type delta = Type(1.e9);
			size_t iteration = 0;
			while (std::abs(delta) > TOL)
			{
				if (iteration++ > MAX_NUM_ITERATIONS) break;

				// Calculate function and it derivative.
				const Type tV = std::tanh(x - yV);
				Type F = s + g * tV - yV;
				Type dF = g * (tV * tV - 1) - 1;

				// Get Newton direction
				Type dir = -F / dF;

				// Line search
				size_t iarm = 0;
				Type lambda = 1;
				delta = lambda * dir;
				Type yN = yV + delta;
				while (std::abs(s + g * std::tanh(x - yN) - yN) > (1 - alpha * lambda) * std::abs(F))
				{
					lambda *= sigma;
					delta = lambda * dir;
					yN = yV + delta;

					if (iarm++ > 8) break;
				}

				yV = yN;
			}

			// Update integrator
			s = 2 * yV - s;

			return yV;
		}

		//==============================================================================
		void update() noexcept
		{
			g = std::tan(c<Type>::pi * cutoff / sampleRate);
		}

		// Parameters of a processor
		Type sampleRate{ 44100. }, cutoff{ 200. }, gain{ 1 };
		Type g{};

		std::vector<Type> st{ 2 };

		// Newton parameters
		static constexpr Type alpha = Type(1.e-4);
		static constexpr Type sigma = Type(0.1);
		static constexpr Type TOL = Type(1.e-5);
		static constexpr size_t MAX_NUM_ITERATIONS = 32;
	};
}
