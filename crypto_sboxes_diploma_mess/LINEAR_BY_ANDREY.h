#include <cstdint>
#include <vector>
#include <algorithm>
#include <numeric>
#include <omp.h>
#pragma once
// В зависимости от компилятора использует соответствующую функцию подсчета числа единичных бит.
#if defined(__GNUC__)
#define popcount64 __builtin_popcountll
#define popcount32 __builtin_popcount
#define popcount16 __builtin_popcount
#elif defined(_MSC_VER)
#include <intrin.h>
#define popcount64 __popcnt64
#define popcount32 __popcnt
#define popcount16 __popcnt16
#endif

// Проверка разрядности системы для Visual C/C++.
#if _WIN32 || _WIN64
#if _WIN64
#define _X64
#else
#define _X86
#endif
#endif
// Проверка разрядности системы для GCC (или Clang).
#if __GNUC__
#if __x86_64__ || __ppc64__
#define _X64
#else
#define _X86
#endif
#endif

// Определение размера машинного слова в зависимости от рязрядности системы.
#if defined(_X64)
typedef uint64_t word; // Машинное слово.
#else
typedef uint32_t word;
#endif

typedef uint8_t byte;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;

inline u16 popcount(u16 x) { return popcount16(x); }
inline u32 popcount(u32 x) { return popcount32(x); }
inline u64 popcount(u64 x) { return popcount64(x); }

typedef u16 sub_vector; // Вектор, к которому применяется перестановка.

const int BITS_IN_BYTE = 8;
const int BITS_IN_WORD = BITS_IN_BYTE * sizeof(word);
const word BITS_IN_SUB_VECTOR = sizeof(sub_vector) * BITS_IN_BYTE; // Число битов в векторе.
const word SUB_VECTOR_AMOUNT = 1 << BITS_IN_SUB_VECTOR; // Число всевозможных векторов.
const word SUB_VECTOR_AMOUNT_DIV_BYTE = SUB_VECTOR_AMOUNT / BITS_IN_BYTE; // Число всевозможных векторов, поделённое на длину байта (т.е. на 8).
const word SUB_VECTOR_AMOUNT_DIV_WORD = SUB_VECTOR_AMOUNT_DIV_BYTE / sizeof(word); // Число всевозможных векторов, поделённое на длину машинного слова.
const size_t DOT_MATRIX_SIZE = SUB_VECTOR_AMOUNT * SUB_VECTOR_AMOUNT_DIV_BYTE; // Размер матриц AX и BY (в байтах).
typedef word(*WordDotMatrix)[SUB_VECTOR_AMOUNT_DIV_WORD];

// Скалярное произведение.
inline sub_vector dot(sub_vector x, sub_vector y)
{
	return popcount((sub_vector)(x & y)) & 1;
}

/* Создает двумерный массив (монолитный буфер памяти), содержащий скалярные произведения чисел от 0 до (N - 1) на
перестановленные числа, причем всё это в упакованном виде. */
WordDotMatrix generate_dot_matrix(const sub_vector * substitution, bool invert_matrix)
{
	WordDotMatrix result;
	result = static_cast<decltype(result)>(operator new(DOT_MATRIX_SIZE)); // Выделение памяти для содержания матрицы.

	for (word a = 0; a < SUB_VECTOR_AMOUNT; ++a)
	{
		for (word x = 0; x < SUB_VECTOR_AMOUNT; ++x)
		{
			auto position_in_line = x / BITS_IN_WORD;

			// Обнуление памяти.
			if (x % BITS_IN_WORD == 0)
				result[a][position_in_line] = invert_matrix ? ~(word)0 : 0;

			auto y = substitution[x];
			word bit_value = dot(a, y);
			result[a][position_in_line] ^= (bit_value << (x % BITS_IN_WORD));
		}
	}
	return result;
}

std::vector<sub_vector> get_identical_substitution()
{
	std::vector<sub_vector> result(SUB_VECTOR_AMOUNT);
	std::iota(std::begin(result), std::end(result), 0);
	return result;
}

std::vector<sub_vector> get_random_substitution()
{
	std::vector<sub_vector> result(SUB_VECTOR_AMOUNT);
	std::iota(std::begin(result), std::end(result), 0);
	std::random_shuffle(result.begin(), result.end());
	return result;
}


auto find_linear_characteristic_for_a_b_fast(const word * inverse_AX_matrix_line, const word * BY_matrix_line)
{
	word linear_characteristic_for_specified_a_b = 0;
	for (word i = 0; i < SUB_VECTOR_AMOUNT_DIV_WORD; ++i)
	{
		auto temp = inverse_AX_matrix_line[i] ^ BY_matrix_line[i];
		linear_characteristic_for_specified_a_b += popcount(temp);
	}
	return linear_characteristic_for_specified_a_b;
}
// Находит максимальную линейную характеристику для конкретного a и всех b.
auto find_max_linear_characteristic_for_a(const WordDotMatrix inverse_AX_matrix, const WordDotMatrix BY_matrix, sub_vector a)
{
	const word * inverse_AX_matrix_line = inverse_AX_matrix[a];

	word max_linear_characteristic_for_a = 0;
	for (word b = a != 0 ? 0 : 1; b < SUB_VECTOR_AMOUNT; ++b)
	{
		word linear_characteristic_for_specified_a_b = find_linear_characteristic_for_a_b_fast(inverse_AX_matrix_line, BY_matrix[b]);

		max_linear_characteristic_for_a = std::max(max_linear_characteristic_for_a, linear_characteristic_for_specified_a_b);
		++b;
	}
	return max_linear_characteristic_for_a;
}

extern "C"
size_t find_max_linear_characteristic(const sub_vector * substitution_table)
{
	WordDotMatrix inverse_AX_matrix = generate_dot_matrix(get_identical_substitution().data(), true);
	WordDotMatrix BY_matrix = generate_dot_matrix(substitution_table, false);

	word max_linear_characteristic = 0;
	omp_set_num_threads(7);
#pragma omp parallel for
	for (int a = 0; a < SUB_VECTOR_AMOUNT; ++a)
	{
		auto max_linear_characteristic_for_a = find_max_linear_characteristic_for_a(inverse_AX_matrix, BY_matrix, a);
#pragma omp critical
		{
			max_linear_characteristic = std::max(max_linear_characteristic, max_linear_characteristic_for_a);
		}
	}

	operator delete(inverse_AX_matrix);
	operator delete(BY_matrix);

	return (1 << 15) - max_linear_characteristic;
}