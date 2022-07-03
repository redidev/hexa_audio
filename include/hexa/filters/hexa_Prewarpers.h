#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <type_traits>

namespace hexa
{
	/**
	 * Realization of a classical prewarper for BLT (aka tangent prewarper)
	 */
	template <typename Type, typename = std::enable_if_t<std::is_floating_point_v<std::remove_cv_t<Type>>>>
	class SimplePrewarper
	{
	public:
		SimplePrewarper() = default;

		void setup(Type newSampleRate) noexcept
		{
			assert(newSampleRate > Type(0));
			sampleRate = newSampleRate;
		}

		Type g(Type freq) const noexcept
		{
			assert(freq > 0 && freq < sampleRate / 2);
			return std::tan(c<Type>::pi * freq / sampleRate);
		}

		Type mu(Type freq) const noexcept { return g(freq) * 2 * sampleRate; }

	private:
		Type sampleRate{ 44100. };
	};

	/**
	 * Realization of a smoothed prewarper with a transition point (see Vadim Zavalishin's book)
	 */
	template <typename Type, typename = std::enable_if_t<std::is_floating_point_v<std::remove_cv_t<Type>>>>
	class TaylorPrewarper
	{
	public:
		TaylorPrewarper() noexcept
		{ update(); }

		void setup(Type newSampleRate) noexcept
		{
			assert(newSampleRate > Type(0));
			sampleRate = newSampleRate;
			update();
		}

		void setTransitionFrequency(Type transFreq) noexcept
		{
			transPoint = std::clamp(transFreq, Type(16.e3), Type(20.e3));
			update();
		}

		Type g(Type freq) const noexcept
		{
			return freq < transPoint ? std::tan(c<Type>::pi * freq / sampleRate) : a * c<Type>::twoPi * freq + b;
		}

		Type mu(Type freq) const noexcept { return  g(freq) * 2 * sampleRate; }

	private:
		void update() noexcept
		{
			Type A = std::tan(c<Type>::pi * transPoint / sampleRate);
			a = (1 + A * A) / (2 * sampleRate);
			b = A - a * c<Type>::twoPi * transPoint;
		}

		Type sampleRate{ 44100. }, transPoint{ 16.5e3 };
		Type a{}, b{};
	};
}