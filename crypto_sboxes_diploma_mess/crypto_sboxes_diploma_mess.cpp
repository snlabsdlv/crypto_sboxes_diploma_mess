// veligura.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <fstream>
#include <memory>
#include <bitset>
#include <omp.h>
#include <algorithm>
#include <stdio.h>
#include <tchar.h>

#include "Kuzn.h"
#include "AES.h"
#include "AES16x16.h"
#include "Butterfly.h"
#include "StandartSboxes.h"
#include <intrin.h>
#include <vector>
#include <algorithm>
#include <bitset>
//#include "LINEAR_BY_ANDREY.h"
#define SIZE 8
using namespace std;

//Ќ≈ѕ–ј¬»Ћ№Ќќ
void avalanche(int *sub) {
	int counter = 0;
	float sum = 0;
	for (int i = 0; i < SIZE; i++) {
		for (int x = 0; x < 1 << SIZE; x++) {
			sum += bitset<SIZE>(sub[x] ^ sub[x ^ (1 << i)]).count();
		}
	}
	cout << sum / ((1 << SIZE)*SIZE) << endl;
}

void SAC(int *sub) {
	int sums[SIZE] = { 0 };
	uint64_t sum = 0;
	for (int i = 0; i < SIZE; i++) {
		fill(begin(sums), begin(sums) + SIZE, 0);
		for (int x = 0; x < 1 << SIZE; x++) {
			for (int j = 0; j < SIZE; j++) {
				sums[j] += bitset<SIZE>(sub[x] ^ sub[x ^ (1 << i)])[j];
			}
		}
		//cout << "i = " << i << endl;
		cout << "(";
		for (int j = 0; j < SIZE; j++) {
			//j!=(SIZE-1) ? cout << sums[j] << ", " : cout << sums[j] << ")";
			sum += sums[j];
		}
	}
	cout << sum / (SIZE * SIZE) << endl;
}
void SAC2(int *sub)
{
	for (int i = 0; i < SIZE; i++)
	{

		int bin[SIZE] = { 0 };
		int bit = 1 << i;
		printf("bit : %d\n", i);
		int z = 0;
		int count = 0;
		for (int x = 0; x < (1 << SIZE); x++)
		{
			int a = sub[x];
			//printf("a: %X\n", a); 
			int b = sub[x^bit];
			// printf("b: %X\n", b); 
			int c = a ^ b;
			// printf("c: %X\n", c); 
			for (int j = 0; j < SIZE; j++)
			{
				int d = ((1 << j) & c) >> j;
				bin[j] += d;
			}
			/* for (int k = 0; k < SIZE; k++)
			{
			printf("%d\n", bin[k]);
			}
			_getch();
			*/

			//count++; 

		}
		// printf("count: %d\n", count); 
		for (int k = 0; k < SIZE; k++)
		{
			printf("%d\n", bin[k]);
		}
	}
}
int* create_set_of_zhegalkin_pol(int* sub) {
	int* pols = (int*)malloc((1 << SIZE) * sizeof(int));
	int* table_copy = (int*)malloc((1 << SIZE) * sizeof(int));

	memcpy(table_copy, sub, (1 << SIZE) * sizeof(int));
	pols[0] = table_copy[0];

	for (int i = 1; i < (1 << SIZE); i++)
	{
		int buf[1 << SIZE];
		memcpy(buf, table_copy, (1 << SIZE) * sizeof(int));

		for (int j = ((1 << SIZE) - i); j > 0; j--)
		{
			buf[j - 1] = table_copy[j] ^ table_copy[j - 1];
		}
		pols[i] = buf[0];
		memcpy(table_copy, buf, (1 << SIZE) * sizeof(int));
	}
	return pols;
}

int* calc_weights(int* sub) {
	int *weights = (int*)malloc(SIZE * sizeof(int));
	for (int i = 0; i < SIZE; i++)
		weights[i] = 0;

	for (int i = 0; i < (1 << SIZE); i++)
		for (int j = 0; j < SIZE; j++)
			weights[j] += (sub[i] & (1 << j)) >> j;
	return weights;
}

int* calc_monoms(int* zhegals) {
	int *nums_of_monoms = (int*)malloc(SIZE * sizeof(int));
	for (int i = 0; i < SIZE; i++)
		nums_of_monoms[i] = 0;
	for (int i = 0; i < (1 << SIZE); i++)
		for (int j = 0; j < SIZE; j++)
			nums_of_monoms[j] += (zhegals[i] & (1 << j)) >> j;
	return nums_of_monoms;
}

int calc_all_monoms(int* zhegals) {
	int num_of_monoms = 0;
	int pol[1 << SIZE] = { 0 };
	for (int i = 0; i < (1 << SIZE); i++)
		for (int j = 0; j < SIZE; j++)
			pol[i] |= (zhegals[i] & (1 << j)) >> j;
	for (int i = 0; i < (1 << SIZE); i++)
		num_of_monoms += pol[i];
	return num_of_monoms;
}

int* calc_fourier_spectrum(int* sub) {
	//0-127 спектр 1-й координатной функции,128-255 - 2-й и т.д. 
	int *Cf = (int*)malloc(SIZE * (1 << SIZE) * sizeof(int));
	for (int i = 0; i < SIZE * (1 << SIZE); i++)
		Cf[i] = 0;
	for (int i = 0; i < SIZE; i++)
	{
		for (int a = 0; a < (1 << SIZE); a++)
		{
			for (int x = 0; x < (1 << SIZE); x++)
			{
				//a1x1+a2x2+...+anxn проверка четности суммы
				int sum = 0;
				for (int i = 0; i < SIZE; i++)
					sum += ((a&x) >> i) & 1;
				if ((sum & 1) == 0)
					Cf[a + i * (1 << SIZE)] += (sub[x] & (1 << i)) >> i;
				else
					Cf[a + i * (1 << SIZE)] -= (sub[x] & (1 << i)) >> i;
			}
		}
	}
	return Cf;
}
//Cf-спектр фурье
int* calc_walsh_spectrum(int *Cf)
{
	int *Zf = (int*)malloc(SIZE * (1 << SIZE) * sizeof(int));
	for (int i = 0; i < SIZE; i++)
		for (int a = 0; a < (1 << SIZE); a++)
		{
			if (a == 0)
				Zf[a + i * (1 << SIZE)] = (1 << SIZE) - 2 * Cf[a + i * (1 << SIZE)];
			else
				Zf[a + i * (1 << SIZE)] = -2 * Cf[a + i * (1 << SIZE)];
		}
	return Zf;
}
//zhegalkin-массив полиномов ∆егалкина координатных функций
int non_linearity_degree(int* zhegalkin, int pos)
{
	int maxDegree = 0;

	for (int i = 0; i < (1 << SIZE); i++)
	{
		if ((zhegalkin[i] & (1 << pos)) != 0)
		{
			int curDegree = 0;
			for (int j = 0; j < SIZE; j++)
			{
				if ((i & (1 << j)) != 0)
					curDegree++;
			}
			if (curDegree > maxDegree)
				maxDegree = curDegree;
		}
	}

	return maxDegree;
}
int numerate_composition(int *sub)//sub-массив жегалкинов
{
	int minDegree = SIZE, degree;
	for (int i = 1; i < (1 << SIZE); i++)
	{
		int compositon[1 << SIZE] = { 0 };
		for (int j = 0; j < SIZE; j++)
		{
			if ((i & 1 << j) != 0)
				for (int k = 0; k < 1 << SIZE; k++)
					compositon[k] ^= (sub[k] & (1 << j)) >> j;
		}
		degree = non_linearity_degree(compositon, 0);
		{
			if (minDegree > degree)
				minDegree = degree;
		}
	}
	return minDegree;
}
/*
int absolute_bias(int* sub)
{
	int lbgS[(1 << SIZE)][1 << SIZE];
	int laS[(1 << SIZE)][1 << SIZE];
	for (int i = 1; i < (1 << SIZE); i++)
	{
		for (int k = 0; k < (1 << SIZE); k++)
		{
			int la = 0;
			int lbg = 0;
			int consta = (i & k);
			int constb = (i & sub[k]);
			for (int q = 0; q < SIZE; q++)
			{
				la += (consta & (1 << q)) >> q;
				lbg += (constb & (1 << q)) >> q;
			}
			laS[i][k] = la % 2;
			lbgS[i][k] = lbg % 2;
		}
	}

	int max_bias = -1;
	for (int i = 1; i < (1 << SIZE); i++)
	{
		for (int j = 1; j < (1 << SIZE); j++)
		{
			int counter = 0;
			for (int k = 0; k< (1 << SIZE); k++)
			{
				if (laS[i][k] == lbgS[j][k])
					counter++;
			}
			int corr = counter;
			if (corr > max_bias)
				max_bias = corr;
		}
	}
	return max_bias - (1 << (SIZE-1));
}
*/
int my_bias(int* sub) {
	int maxValue = 0;
	int value = 0;
	for (int i = 1; i < (1 << SIZE); i++) {
		for (int j = 1; j < (1 << SIZE); j++) {
			value = 0;
			for (int x = 0; x < (1 << SIZE); x++) {
				value += !(__popcnt((x & i) ^ (sub[x] & j)) & 1);
			}
			maxValue = maxValue < value ? value : maxValue;
		}
	}
	return maxValue - (1 << (SIZE - 1));
}

int alfa_beta_test(int* sub)
{
	int maxValue = 0;
	//uint64_t maxSum = 0;
	//uint64_t sum = 0;
	//uint64_t sums[1 << SIZE] = { 0 };
	for (int a = 1; a < (1 << SIZE); a++)
	{
		int b[1 << SIZE] = { 0 };
		for (int x = 0; x < (1 << SIZE); x++)
		{
			if ((b[sub[x] ^ sub[x^a]]++) > maxValue)
				maxValue = b[sub[x] ^ sub[x^a]];
		}

		//sum = 0;
		//for (int k = 0; k < (1 << SIZE); k++) {
		//	 sum += pow(b[k], 9);
		//}
		//maxSum = std::max(sum, maxSum);

	}
	/*
	for (int k = 0; k < (1 << SIZE); k++) {
		maxSum = maxSum < sums[k] ? sums[k] : maxSum;
	}
	*/
	//printf("%lu \n", maxSum);
	return maxValue;
}

int* calc_fourier_spectrum1(int* sub) {
	int *Cf = (int*)malloc((1 << SIZE) * sizeof(int));
	for (int i = 0; i < (1 << SIZE); i++)
		Cf[i] = 0;
	for (int a = 0; a < (1 << SIZE); a++)
	{
		for (int x = 0; x < (1 << SIZE); x++)
		{
			int sum = 0;
			for (int i = 0; i < SIZE; i++)
				sum += ((a&x) >> i) & 1;
			if ((sum & 1) == 0)
				Cf[a] += sub[x] & 1;
			else
				Cf[a] -= sub[x] & 1;
		}
	}

	return Cf;
}
int* calc_walsh_spectrum1(int *Cf)
{
	int *Zf = (int*)malloc((1 << SIZE) * sizeof(int));
	for (int a = 0; a < (1 << SIZE); a++)
	{
		if (a == 0)
			Zf[a] = (1 << SIZE) - 2 * Cf[a];
		else
			Zf[a] = -2 * Cf[a];
	}
	return Zf;
}

void check_perf_K() {
	uint8_t plain[16] = { 0 };
	uint8_t plain_xor_1[16] = { 0 };

	uint8_t out[16] = { 0 };
	uint8_t out_xor_1[16] = { 0 };
	uint8_t keys[16] = { 0 };

	bool isPerf[128] = { false };
	bool allIsPerf[128] = { false };
	bool perfect = false;
	for (int blocki = 0; blocki < 16; blocki++)
	{
		for (int coordi = 0; coordi < 8; coordi++)
		{
			for (int m = 0; m < 128; m++) {
				isPerf[m] = false;
			}
			for (int i = 0; i < 256; i++) {
				for (int k = 0; k < 16; k++) {
					plain[k] = 0;
					plain_xor_1[k] = 0;
				}
				plain[blocki] = i;
				plain_xor_1[blocki] = i ^ (1 << coordi);

				encryptBlockWithGost15(keys, plain);
				encryptBlockWithGost15(keys, plain_xor_1);

				for (int j = 0; j < 16; j++) {
					for (int k = 0; k < 8; k++) {
						if ((plain[j] & (1 << k)) != (plain_xor_1[j] & (1 << k))) {
							isPerf[j * 8 + k] = true;
						}
					}
				}
				allIsPerf[coordi + blocki * 8] = true;
				for (int m = 0; m < 128; m++) {
					if (isPerf[m] == false) {
						allIsPerf[coordi + blocki * 8] = false;
					}
				}
				if (allIsPerf[coordi + blocki * 8]) {
					break;
				}
			}
		}
	}
	perfect = true;
	for (int m = 0; m < 128; m++) {
		if (allIsPerf[m] == false) {
			perfect = false;
		}
	}
	if (perfect) {
		cout << "PERFECT\n";
	}
}

void check_perf_AES() {
	uint8_t plain[16] = { 0 };
	uint8_t plain_xor_1[16] = { 0 };

	uint8_t out[16] = { 0 };
	uint8_t out_xor_1[16] = { 0 };
	uint8_t keys[16] = { 0 };

	bool isPerf[128] = { false };
	bool allIsPerf[128] = { false };
	bool perfect = false;
	for (int blocki = 0; blocki < 16; blocki++)
	{
		for (int coordi = 0; coordi < 8; coordi++)
		{
			for (int m = 0; m < 128; m++) {
				isPerf[m] = false;
			}
			for (int i = 0; i < 256; i++) {
				plain[blocki] = i;
				plain_xor_1[blocki] = i ^ (1 << coordi);

				AESEncrypt(plain, keys, out);
				AESEncrypt(plain_xor_1, keys, out_xor_1);

				for (int j = 0; j < 16; j++) {
					for (int k = 0; k < 8; k++) {
						if ((out[j] & (1 << k)) != (out_xor_1[j] & (1 << k))) {
							isPerf[j * 8 + k] = true;
						}
					}
				}
				allIsPerf[coordi + blocki * 8] = true;
				for (int m = 0; m < 128; m++) {
					if (isPerf[m] == false) {
						allIsPerf[coordi + blocki * 8] = false;
					}
				}
				if (allIsPerf[coordi + blocki * 8]) {
					break;
				}
			}
		}
	}
	perfect = true;
	for (int m = 0; m < 128; m++) {
		if (allIsPerf[m] == false) {
			perfect = false;
		}
	}
	if (perfect) {
		cout << "PERFECT\n";
	}
}

#define ROTL8(x,shift) ((uint8_t) ((x) << (shift)) | ((x) >> (8 - (shift))))
void initialize_aes_sbox(uint8_t sbox[256]) {
	uint8_t p = 1, q = 1;

	/* loop invariant: p * q == 1 in the Galois field */
	do {
		/* multiply p by 3 */
		p = p ^ (p << 1) ^ (p & 0x80 ? 0x1B : 0);

		/* divide q by 3 (equals multiplication by 0xf6) */
		q ^= q << 1;
		q ^= q << 2;
		q ^= q << 4;
		q ^= q & 0x80 ? 0x09 : 0;

		/* compute the affine transformation */
		uint8_t xformed = q ^ ROTL8(q, 1) ^ ROTL8(q, 2) ^ ROTL8(q, 3) ^ ROTL8(q, 4);

		sbox[p] = xformed ^ 0x63;
	} while (p != 1);

	/* 0 is a special case since it has no inverse */
	sbox[0] = 0x63;
}


int main()
{
	/*
	uint8_t qwe[256] = { 0 };
	initialize_aes_sbox(qwe);
	for (uint8_t q = 0; q < 256; q++) {
		printf("%d %d \n", qwe[q], q ? q ^ ROTL8(q, 1) ^ ROTL8(q, 2) ^ ROTL8(q, 3) ^ ROTL8(q, 4) ^ 0x63 : 0x63);
	}
	*/

	int gsbox[1 << SIZE] = { 0 };
	for (unsigned int i = 0; i < 1 << SIZE; i++) {
		gsbox[i] = speck2(i);
	}


	/*
	for (int j = 1; j < 256 - 2; j++) {
		uint16_t zz[256] = { 0 };
		for (int in = 0; in < 256; in++) {
			zz[in] = in;
			for (int i = 0; i < j; i++) {
				zz[in] = multiplyInGF256_0x1b(zz[in], in);
			}
		}
		printf("uint16_t g%d_256[256] = {",j+1);
		for (int a = 0; a < 256; a++) {
			printf("%d, ", zz[a]);
		}
		cout << "};\n";
	}
	*/


	uint16_t x, inv;
	int S[0x10000];
	//uint16_t SINV[0x10000];
	int i, j;
	// Compute the S-box values
	for (i = 0; i <= 0xFFFF; i++)
	{
		x = (uint16_t)i;
		S[i] = forward(x);
		//SINV[i] = inverse(x);
		//if (inverse(forward(x)) != x)
		//{
		//	printf("FAILURE with 0x%x\n", x);
		//	return -1;
		//}
	}

	//printf("MAXVALUE=%d\n",alfa_beta_test(gsbox));
	//alfa_beta_test(gsbox);

	//////////////// int butt8table[1 << SIZE] = { 0 };
	/*
	for (int pi2index = 0; pi2index < 13; pi2index++) {
		for (int pi1index = 0; pi1index < 13; pi1index++) {
			int butt8table[1 << SIZE] = { 0 };
			for (int i = 0; i < 1 << SIZE; i++) {
				butt8table[i] = butt8Jimenez(i, gs[pi1index], gs[pi2index]);
			}
			printf("pi1 = %d, pi2 = %d: %d\n", pi1index + 2, pi2index + 2, alfa_beta_test(butt8table));
		}
	}
	*/
	/*
	for (int i = 0; i < 1 << SIZE; i++) {
		butt8table[i] = butt8FominB2(i);
	}
	printf("%d\n", alfa_beta_test(butt8table));
	*/
	/*
	int butt16table[1 << SIZE] = { 0 };
	ofstream wfile;
	int *polynomials;
	int polDegree;
	omp_set_num_threads(7);
#pragma omp parallel for private(butt16table, polynomials, polDegree)
	for (int pi1index = 0; pi1index < 253; pi1index++) {
		for (int i = 0; i < 1 << SIZE; i++) {
			butt16table[i] = butt16Fomin(i, g256s[pi1index]);
		}
		int ps = alfa_beta_test(butt16table);
		if (ps <= 20) {
			int minDegree = 999;
			polynomials = create_set_of_zhegalkin_pol(butt16table);
			for (int i = 0; i < SIZE; i++)
			{
				polDegree = non_linearity_degree(polynomials, i);
				if (polDegree < minDegree) {
					minDegree = polDegree;
				}
			}

			printf("pi1 = %d, pi2 = 254, nonlin = %d, Ps = %d\n", pi1index + 2, minDegree, ps);
#pragma omp critical
			{
				wfile.open("16Fomin_fixPi2(254)_0x1b.txt", ios::app);
				wfile << "pi1 = " << pi1index + 2 << ", pi2 = " << 254 << " nonlin = " << minDegree << " Ps = " << ps << endl;
				wfile.close();
			}
		}
	}
	*/
	/*
	int butt16table[1 << SIZE] = { 0 };
	int *polynomials;
	int polDegree;
	int indexes_plus_2[16] = { 28, 37, 56, 73, 74, 131, 146, 148, 164, 191, 193, 239, 247, 251, 253, 254 };
	omp_set_num_threads(7);
#pragma omp parallel for private(butt16table, polynomials, polDegree)
	for (int pi1index = 0; pi1index < 16; pi1index++) {
		for (int i = 0; i < 1 << SIZE; i++) {
			butt16table[i] = butt16Jimenez(i, g256s[indexes_plus_2[pi1index] - 2], g253_256);
		}
		int minDegree = 999;
		polynomials = create_set_of_zhegalkin_pol(butt16table);
		for (int i = 0; i < SIZE; i++)
		{
			polDegree = non_linearity_degree(polynomials, i);
			if (polDegree < minDegree) {
				minDegree = polDegree;
			}
		}
#pragma omp critical
		{
			printf("g %d = %d\n", indexes_plus_2[pi1index], minDegree);
		}
	}
	*/

	int butt16table[1 << SIZE] = { 0 };
	for (int i = 0; i < 1 << SIZE; i++) {
		butt16table[i] = butt16Jimenez(i, g254_256, g254_256);
		//butt16table[i] = (AES_sbox[(i&0xff00)>>8]<<8) ^ AES_sbox[i&0xff];
		//butt16table[i] = butt16from4_2(i, gs[k1], gs[k2]);
	}
	//printf("maxlin = %d \n", find_max_linear_characteristic(butt16table));
	//printf("maxlin = %d \n", absolute_bias(Skipjack_sbox));
	//cout << endl;
	//printf("%d\n", my_bias(Kuznechik_sbox));
	printf("%d\n", alfa_beta_test(MenyachikhinBeLT));

	int maxW = -1;
	int *spectrum = calc_walsh_spectrum(calc_fourier_spectrum(MenyachikhinBeLT));
	for (int i = 0; i < SIZE*(1 << SIZE); i++) {
		maxW = max(maxW, spectrum[i]);
	}
	int nonlinearity = (1 << (SIZE - 1)) - maxW / 2;
	printf("maxW = %d NONLINEARITY= %d", maxW, nonlinearity);

	//SAC(Skipjack_sbox);
	//avalanche(gsbox);
	return 0;
}

