#pragma once
#include <stdint.h>
#include <string.h>
#include "Kuzn_tables.h"
#include "gSbox.h"
#include "Galois256_65536.h"
// NumberOfRounds = 2 (==1 round)
enum {
	NumberOfRounds = 10,
	NumberOfRoundsInKeySchedule = 8,

	BlockLengthInBytes = 128 / 8,
	KeyLengthInBytes = 256 / 8,
};
const size_t WorkspaceOfScheduleRoundKeys = BlockLengthInBytes * 2;
static uint8_t multiplyInGF256(
	uint8_t left,
	uint8_t right
) {
	if (left && right) {
		int productLogarithm_ = 0;
		uint8_t leftLogarithm_ = logarithmicTable[left];
		uint8_t rightLogarithm_ = logarithmicTable[right];
		productLogarithm_ = leftLogarithm_ + rightLogarithm_;
		if (productLogarithm_ > 0xff) { productLogarithm_ -= 0xff; }
		return exponentialTable[productLogarithm_];
	}
	else {
		return 0;
	}
}

static void applyXTransformation(
	const uint8_t * key,
	const uint8_t *input,
	uint8_t * output
) {
	uint64_t key_[2], input_[2], output_[2];

	/* Compiler usually sees this through. Advantages of using memcpy:  */
	memcpy(key_, key, BlockLengthInBytes);
	memcpy(input_, input, BlockLengthInBytes);

	output_[0] = input_[0] ^ key_[0];
	output_[1] = input_[1] ^ key_[1];

	memcpy(output, output_, BlockLengthInBytes);
}

static void applySTransformation(
	uint8_t *block
) {
	for (int byteIndex_ = 0; byteIndex_ < BlockLengthInBytes; byteIndex_ += 8) {
		block[byteIndex_ + 0] = Pi[block[byteIndex_ + 0]];
		block[byteIndex_ + 1] = Pi[block[byteIndex_ + 1]];
		block[byteIndex_ + 2] = Pi[block[byteIndex_ + 2]];
		block[byteIndex_ + 3] = Pi[block[byteIndex_ + 3]];

		block[byteIndex_ + 4] = Pi[block[byteIndex_ + 4]];
		block[byteIndex_ + 5] = Pi[block[byteIndex_ + 5]];
		block[byteIndex_ + 6] = Pi[block[byteIndex_ + 6]];
		block[byteIndex_ + 7] = Pi[block[byteIndex_ + 7]];
	}
}

static void applySTransformation16(
	unsigned char *block
) {
	for (int byteIndex_ = 0; byteIndex_ < 16; byteIndex_ += 2) {
		speck2ptr(block + byteIndex_);
	}
}

static void applyInversedSTransformation(
	uint8_t *block
) {
	for (int byteIndex_ = 0; byteIndex_ < BlockLengthInBytes; byteIndex_ += 8) {
		block[byteIndex_ + 0] = InversedPi[block[byteIndex_ + 0]];
		block[byteIndex_ + 1] = InversedPi[block[byteIndex_ + 1]];
		block[byteIndex_ + 2] = InversedPi[block[byteIndex_ + 2]];
		block[byteIndex_ + 3] = InversedPi[block[byteIndex_ + 3]];

		block[byteIndex_ + 4] = InversedPi[block[byteIndex_ + 4]];
		block[byteIndex_ + 5] = InversedPi[block[byteIndex_ + 5]];
		block[byteIndex_ + 6] = InversedPi[block[byteIndex_ + 6]];
		block[byteIndex_ + 7] = InversedPi[block[byteIndex_ + 7]];
	}
}


static void applyLTransformation(
	const uint8_t *input,
	uint8_t *output
) {
	for (int byteIndex_ = 0; byteIndex_ < BlockLengthInBytes; ++byteIndex_) {
		uint8_t cache_ = 0;
		for (int addendIndex_ = 0; addendIndex_ < BlockLengthInBytes; ++addendIndex_) {
			cache_ ^= multiplyInGF256(LTransformationMatrix[addendIndex_][byteIndex_],
				input[addendIndex_]);
		}
		output[byteIndex_] = cache_;
	}
}

static void applyLTransformation16(
	const uint16_t *input,
	uint16_t *output
) {
	for (int byteIndex_ = 0; byteIndex_ < 8; ++byteIndex_) {
		uint16_t cache_ = 0;
		for (int addendIndex_ = 0; addendIndex_ < 8; ++addendIndex_) {
			cache_ ^= multiplyInGF65536(LTransformationMatrix16[addendIndex_][byteIndex_],
				input[addendIndex_]);
		}
		output[byteIndex_] = cache_;
	}
}

static void applyInversedLTransformation(
	const uint8_t *input,
	uint8_t *output
) {
	for (int byteIndex_ = 0; byteIndex_ < BlockLengthInBytes; ++byteIndex_) {
		uint8_t cache_ = 0;
		for (int addendIndex_ = 0; addendIndex_ < BlockLengthInBytes; ++addendIndex_) {
			cache_ ^= multiplyInGF256(inversedLTransformationMatrix[addendIndex_][byteIndex_],
				input[addendIndex_]);
		}
		output[byteIndex_] = cache_;
	}
}


static void applyXSLTransformation(
	const uint8_t *key,
	uint8_t *block,
	uint8_t *temporary
) {
	applyXTransformation(key, block, temporary);
	//applySTransformation(temporary);
	applySTransformation16(temporary);
	//applyLTransformation(temporary, block);
	applyLTransformation16((uint16_t *)temporary, (uint16_t *)block);
}


static void applyInversedSLXTransformation(
	const uint8_t *key,
	uint8_t *block,
	uint8_t *temporary
) {
	applyXTransformation(key, block, temporary);
	applyInversedLTransformation(temporary, block);
	applyInversedSTransformation(block);
}


static void swapBlocks(
	uint8_t * left,
	uint8_t * right,
	uint8_t * temporary
) {
	/* left != right != temp shall hold. */
	memcpy(temporary, left, BlockLengthInBytes);
	memcpy(left, right, BlockLengthInBytes);
	memcpy(right, temporary, BlockLengthInBytes);
}


static void applyFTransformation(
	const uint8_t * key,
	uint8_t * left,
	uint8_t * right,
	uint8_t * temporary1,
	uint8_t * temporary2
) {
	memcpy(temporary1, left, BlockLengthInBytes);
	applyXSLTransformation(key, temporary1, temporary2);
	applyXTransformation(temporary1, right, right);
	swapBlocks(left, right, temporary2);
}


static void fetchKeyScheduleRoundConstant(
	uint8_t index,
	uint8_t * roundConstant,
	uint8_t * temporary
) {
	memset(temporary, 0, BlockLengthInBytes);
	temporary[BlockLengthInBytes - 1] = index;
	applyLTransformation(temporary, roundConstant);
}


static void scheduleRoundKeys(
	uint8_t * roundKeys,
	const void * key,
	uint8_t * memory
) {
	uint8_t *memory_ = memory, *roundKeys_ = roundKeys;

	memcpy(&roundKeys_[0], key, BlockLengthInBytes * 2);

	for (int nextKeyIndex_ = 2, constantIndex_ = 0;
		nextKeyIndex_ != NumberOfRounds;
		nextKeyIndex_ += 2) {
		memcpy(&roundKeys_[BlockLengthInBytes * (nextKeyIndex_)],
			&roundKeys_[BlockLengthInBytes * (nextKeyIndex_ - 2)],
			BlockLengthInBytes * 2);

		for (int feistelRoundIndex_ = 0; feistelRoundIndex_ < NumberOfRoundsInKeySchedule; ++feistelRoundIndex_) {
			applyFTransformation(&roundConstants[BlockLengthInBytes * constantIndex_++],
				&roundKeys_[BlockLengthInBytes * (nextKeyIndex_)],
				&roundKeys_[BlockLengthInBytes * (nextKeyIndex_ + 1)],
				&memory_[0],
				&memory_[BlockLengthInBytes]);
		}
	}
}


void scheduleEncryptionRoundKeysForGost15(
	uint8_t * roundKeys,
	const void * key,
	uint8_t * memory
) {
	scheduleRoundKeys(roundKeys, key, memory);
}


void scheduleDecryptionRoundKeysForGost15(
	uint8_t * roundKeys,
	const void * key,
	uint8_t * memory
) {
	scheduleRoundKeys(roundKeys, key, memory);
}


void encryptBlockWithGost15(
	const uint8_t * roundKeys,
	uint8_t * block
) {
	const uint8_t *roundKeys_ = roundKeys;
	uint8_t cache_[BlockLengthInBytes] = { 0 };
	int round_ = 0;
	//for (; round_ < NumberOfRounds - 1; ++round_) {
	for (; round_ < NumberOfRounds - 1; ++round_) {
		applyXSLTransformation(&roundKeys_[BlockLengthInBytes * round_], block, cache_);
	}
	applyXTransformation(&roundKeys_[BlockLengthInBytes * round_], block, block);
}


void decryptBlockWithGost15(
	const uint8_t * roundKeys,
	uint8_t * block
) {
	const uint8_t *roundKeys_ = roundKeys;
	uint8_t cache_[BlockLengthInBytes] = { 0 };
	int round_ = NumberOfRounds - 1;

	for (; round_ > 0; --round_) {
		applyInversedSLXTransformation(&roundKeys_[BlockLengthInBytes * round_], block, cache_);
	}
	applyXTransformation(&roundKeys_[BlockLengthInBytes * round_], block, block);
}