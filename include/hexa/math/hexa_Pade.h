#pragma once

namespace hexa::pade
{
	//==============================================================================
	template <typename T>
	constexpr T tan(const T& x) noexcept
	{
		auto x2 = x * x;
		return x * (T(945) + (x2 - T(105)) * x2) / (T(945) + (T(15) * x2 - T(420)) * x2);
	}

	//==============================================================================
	template <typename T>
	constexpr T tan2(const T& x) noexcept
	{
		auto x2 = x * x;
		return (T(105) - T(10) * x2) * x / (T(105) + (x2 - T(45)) * x2);
	}

	//==============================================================================
	template <typename T>
	constexpr T tanh(const T& x) noexcept
	{
		auto x2 = x * x;
		return (T(945) + (x2 + T(105)) * x2) * x / (T(945) + (T(15) * x2 + T(420)) * x2);
	}

	template <typename T>
	constexpr T tanh2(const T& x) noexcept
	{
		auto x2 = x * x;
		return (T(10) * x2 + T(105)) * x / (T(105) + (x2 + T(45)) * x2);
	}

	//==============================================================================
	template <typename T>
	constexpr T tanhXdX(const T& x) noexcept
	{
		auto x2 = x * x;
		return (T(945) + (x2 + T(105)) * x2) / (T(945) + (T(15) * x2 + T(420)) * x2);
	}

	template <typename T>
	constexpr T tanhXdX2(const T& x) noexcept
	{
		auto x2 = x * x;
		return (T(10) * x2 + T(105)) / (T(105) + (x2 + T(45)) * x2);
	}

	//==============================================================================
	template <typename T>
	constexpr T asinh(const T& x) noexcept
	{
		auto x2 = x * x;
		return (T(2552051040) + (T(504949578) * x2 + T(2596768160)) * x2) * x / (T(2552051040) 
			+ (T(3022110000) + (T(23477725) * x2 + T(817230750)) * x2) * x2);
	}

	//==============================================================================
	template <typename T>
	constexpr T asinhXdX(const T& x) noexcept
	{
		auto x2 = x * x;
		return (T(2552051040) + (T(504949578) * x2 + T(2596768160)) * x2) / (T(2552051040) 
			+ (T(3022110000) + (T(23477725) * x2 + T(817230750)) * x2) * x2);
	}
}