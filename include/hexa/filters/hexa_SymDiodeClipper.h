#pragma once

#include <vector>
#include <algorithm>
#include <cmath>

#include "../math/hexa_Constants.h"

namespace hexa
{
	/** Implementation of a simple symmetrical diode clipped. */
	template <typename Type>
	class SymDiodeClipper
	{
	public:
		SymDiodeClipper() = default;

		//==============================================================================
		void setFrequency(Type freqHz) noexcept
		{
			cutoff = std::clamp(freqHz, Type(20), Type(16.e3));
			update();
		}

		void setGain(Type gainDb) noexcept
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
			assert(nChans <= st.size());

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
		Type tick(const Type& in, size_t ch) noexcept
		{
			auto&& s = st[ch];

			const Type p = G * (in * gain - s) + s;

			// Capped Newton as described in DAFX-2015 paper (see Ben Holmes)
			// Set initial guess, step size and iterations counter.
			Type y = a * std::asinh(p / b);
			Type delta = 1.e6;
			size_t itr = 0;
			while (std::abs(delta) > TOL)
			{
				if (itr++ > MAX_NUM_ITERATIONS) break;

				const Type F = p - b * std::sinh(y * aInv) - y;
				const Type dF = -1 - b * aInv * std::cosh(y * aInv);

				// Capped step
				delta = std::clamp(-F / dF, -deltaLim, deltaLim);

				y += delta;
			}

			// Update integrator
			s = 2 * y - s;

			return y;
		}

		void update() noexcept
		{
			Type g = std::tan(c<Type>::pi * cutoff / sampleRate);
			G = g / (1 + g);

			a = mu * Vt; aInv = 1 / a;
			b = 2 * R * G * Is;

			deltaLim = a * std::acosh(a / b);
		}

		// Parameters of a processor
		Type sampleRate{ 44100. }, cutoff{ 200. }, gain{ 1. };
		Type G{}, a{}, aInv{}, b{}, deltaLim{};

		std::vector<Type> st{};

		// Parameters for a germanium diode
		static constexpr Type Is = Type(2.52e-9);
		static constexpr Type mu = Type(1.752);		// Ideality
		static constexpr Type Vt = Type(26e-3);

		// Resistor parameters
		static constexpr Type R = Type(2.2e3);

		// NR parameters
		static constexpr Type TOL = Type(1. / 65535);
		static constexpr size_t MAX_NUM_ITERATIONS = 16;
	};
}
