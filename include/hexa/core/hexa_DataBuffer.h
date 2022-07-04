#pragma once

#include <algorithm>
#include <cassert>
#include <iterator>
#include <vector>

namespace hexa
{
	template <typename Type, typename Alloc = std::allocator<Type>>
	class DataBuffer
	{
	public:
		DataBuffer(const DataBuffer& other) = delete;
		DataBuffer& operator= (const DataBuffer& other) = delete;

		DataBuffer(DataBuffer&& other) = default;
		DataBuffer& operator= (DataBuffer&& other) = default;

		DataBuffer(size_t numRows = 32, size_t numCols = 2)
		{
			resize(numRows, numCols);
		}

		void resize(size_t newNumRows, size_t newNumCols)
		{
			numRows = newNumRows; numCols = newNumCols;
			rawData.resize(numRows * numCols, Type(0));
		}

		void clear() noexcept
		{
			std::fill(rawData.begin(), rawData.end(), Type(0));
		}

		void clearColumn(size_t c) noexcept
		{
			std::fill_n(&rawData[flat(0, c)], numRows, Type(0));
		}

		const Type* col(size_t c) const { return &rawData[flat(0, c)]; }

		Type* col(size_t c) { return &rawData[flat(0, c)]; }

		Type& operator() (size_t r, size_t c) { return rawData[flat(r, c)]; }

		const Type& operator() (size_t r, size_t c) const { return rawData[flat(r, c)]; }

		void clearRegion(size_t startRow, size_t length)
		{
			assert(startRow + length < numRows);
			for (size_t c = 0; c < numCols; ++c)
			{
				std::fill_n(&rawData[flat(startRow, c)], length, Type(0));
			}
		}

		size_t getNumRows() const noexcept { return numRows; }

		size_t getNumCols() const noexcept { return numCols; }

		Type* data() { return rawData.data(); }

		const Type* data() const { return rawData.data(); }

	private:
		size_t flat(size_t row, size_t col) const noexcept { return row + numRows * col; }

		std::vector<Type, Alloc> rawData{};
		size_t numRows{}, numCols{};
	};
}