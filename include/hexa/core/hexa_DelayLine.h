#pragma once

#include <vector>
#include <cassert>

#include "hexa_DataBuffer.h"
#include "../math/hexa_Interpolators.h"

namespace hexa
{
	template <typename Type, InterpolationType interp = InterpolationType::CatmullRom, typename Alloc = std::allocator<Type>>
	class DelayLine
	{
	public:
		//==============================================================================
		DelayLine() = delete;

		DelayLine(const DelayLine& other) = delete;
		DelayLine& operator= (const DelayLine& other) = delete;

		DelayLine(DelayLine&& other) = default;
		DelayLine& operator= (DelayLine&& other) = default;

		DelayLine(int reqSize, size_t numChannels = 2)
		{
			resize(reqSize, numChannels);
		}

		//==============================================================================		
		void resize(int newReqSize, size_t newNumChannels)
		{
			assert(newReqSize > 0);
			maxSize = utils::nextPowerOfTwo(newReqSize);
			sizeMsk = maxSize - 1;

			buffer.resize(maxSize, newNumChannels);
			pos.resize(newNumChannels, 0);

			clear();
		}

		void clear() noexcept
		{
			buffer.clear();
			std::fill(pos.begin(), pos.end(), Type(0));
		}

		//==============================================================================
		size_t getMaxSize() const noexcept { return maxSize; }

		size_t getNumChannels() const noexcept { return buffer.getNumCols(); }

		size_t getNumSamples() const noexcept { return buffer.getNumRows(); }

		//==============================================================================
		void push(size_t ch, Type value) noexcept
		{
			auto&& lpos = pos[ch];

			// Write AND Shift
			buffer(lpos, ch) = value;
			lpos = (lpos + 1) & sizeMsk;
		}

		Type operator() (size_t ch, size_t del, double frac = 0) const noexcept
		{
			return interpolate(ch, del, frac);
		}

	private:
		//==============================================================================
		Type interpolate(size_t ch, size_t del, double frac) const noexcept
		{
			if constexpr (interp == InterpolationType::Drop)
			{
				const size_t idx1 = (pos[ch] - del) & sizeMsk;
				return buffer(idx1, ch);
			}
			else if constexpr (interp == InterpolationType::Linear)
			{
				const size_t idx1 = (pos[ch] - del) & sizeMsk;
				const size_t idx2 = (idx1 - 1) & sizeMsk;

				const Type val1 = buffer(idx1, ch);
				const Type val2 = buffer(idx2, ch);

				return op(frac, val1, val2);
			}
			else
			{
				const size_t idx1 = (pos[ch] - del) & sizeMsk;
				const size_t idx2 = (idx1 - 1) & sizeMsk;
				const size_t idx3 = (idx2 - 1) & sizeMsk;
				const size_t idx4 = (idx3 - 1) & sizeMsk;

				const Type val1 = buffer(idx1, ch);
				const Type val2 = buffer(idx2, ch);
				const Type val3 = buffer(idx3, ch);
				const Type val4 = buffer(idx4, ch);

				return op(frac, val1, val2, val3, val4);
			}
		}

		//==============================================================================
		size_t maxSize{}, sizeMsk{};
		std::vector<size_t> pos{};
		DataBuffer<Type, Alloc> buffer{};
		Interpolator<Type, interp> op{};
	};
}