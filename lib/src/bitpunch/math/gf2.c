/*
This file is part of BitPunch
Copyright (C) 2013-2015 Frantisek Uhrecky <frantisek.uhrecky[what here]gmail.com>
Copyright (C) 2013-2015 Andrej Gulyas <andrej.guly[what here]gmail.com>
Copyright (C) 2013-2014 Marek Klein  <kleinmrk[what here]gmail.com>
Copyright (C) 2013-2014 Filip Machovec  <filipmachovec[what here]yahoo.com>
Copyright (C) 2013-2014 Jozef Kudlac <jozef[what here]kudlac.sk>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "gf2.h"
#include "perm.h"

#include <stdlib.h>
#include <bitpunch/debugio.h>
#include <bitpunch/prng/prng.h>

#ifdef BPU_CONF_PRINT
/* ==================================== Print functions ==================================== */
void BPU_printBinaryMsb(uint32_t in, int len) {
	if (len > 0) {
		BPU_printBinaryMsb(in >> 1, len - 1);

		fprintf(stderr, "%d", (int) (in & (0x1u)));
	}
}

void BPU_printBinaryMsbLn(uint32_t in, int len) {
	BPU_printBinaryMsb(in, len);
	fprintf(stderr, "\n");
}

void BPU_printBinaryMsb32(uint32_t in) {
	BPU_printBinaryMsb(in, 32);
}

void BPU_printBinaryMsb32Ln(uint32_t in) {
	BPU_printBinaryMsbLn(in, 32);
}

void BPU_printBinaryLsb(uint32_t in, int len) {
	if (len > 0) {
		fprintf(stderr, "%d", (int) (in & (0x1u)));

		BPU_printBinaryLsb(in >> 1, len - 1);
	}
}

void BPU_printBinaryLsbLn(uint32_t in, int len) {
	BPU_printBinaryLsb(in, len);
	fprintf(stderr, "\n");
}

void BPU_printBinaryLsb32(uint32_t in) {
	BPU_printBinaryLsb(in, 32);
}

void BPU_printBinary32LsbLn(uint32_t in) {
	BPU_printBinaryLsbLn(in, 32);
}

void BPU_printGf2Mat(const BPU_T_GF2_Matrix* m) {
	int i, j, bits_to_print;

	fprintf(stderr, "Matrix size: %dx%d\n", m->k, m->n);

	for (i = 0; i < m->k; i++) {
		fprintf(stderr, "%4d: ",i);

		for (j = 0; j <= m->elements_in_row - 1; j++) {
			if (j == m->elements_in_row-1) {
				if (m->n%(m->element_bit_size) != 0) {
					bits_to_print = m->n % m->element_bit_size;
				}
				else {
					bits_to_print = m->element_bit_size;
				}
			}
			else {
				bits_to_print = m->element_bit_size;
			}
			BPU_printBinaryLsb(m->elements[i][j], bits_to_print);
			// fprintf(stderr, " "); // medzera medzi elementami
		}
		fprintf(stderr, "\n");
	}
}

void BPU_printGf2Vec(const BPU_T_GF2_Vector* v) {
	int j, bits_to_print;

	fprintf(stderr, "Vec (%4d): ", v->len);
	for (j = 0; j <= v->elements_in_row - 1; j++) {
		if (j == v->elements_in_row-1) {
			if (v->len % (v->element_bit_size) != 0) {
				bits_to_print = v->len % v->element_bit_size;
			}
			else {
				bits_to_print = v->element_bit_size;
			}
		}
		else {
			bits_to_print = v->element_bit_size;
		}
		BPU_printBinaryLsb(v->elements[j], bits_to_print);
		fprintf(stderr, " ");
	}
	fprintf(stderr, "\n");
}

void BPU_printGf2VecMsb(const BPU_T_GF2_Vector* v) {
	int j, bits_to_print;

	fprintf(stderr, "Vec (%4d): ", v->len);
	for (j = 0; j <= v->elements_in_row - 1; j++) {
		if (j == v->elements_in_row-1) {
			if (v->len % (v->element_bit_size) != 0) {
				bits_to_print = v->len % v->element_bit_size;
			}
			else {
				bits_to_print = v->element_bit_size;
			}
		}
		else {
			bits_to_print = v->element_bit_size;
		}
		BPU_printBinaryMsbLn(v->elements[j], bits_to_print);
		fprintf(stderr, " ");
	}
	fprintf(stderr, "\n");
}

void BPU_printGf2VecOnes(const BPU_T_GF2_Vector *vec) {
	int i;
	for (i = 0; i < vec->len; ++i)
	{
		if (BPU_gf2VecGetBit(vec, i)) {
			fprintf(stderr, "%d ", i);
		}
	}
	fprintf(stderr, "\n");
}

void BPU_printGf2SparsePoly(const BPU_T_GF2_Sparse_Poly *v) {
	int i;

	fprintf(stderr, "Sparse poly (%i): ", v->weight);
	for (i = 0; i < v->weight; i++) {
		fprintf(stderr, "%3i ", v->index[i]);
	}
	fprintf(stderr, "\n");
}

void BPU_printGf2PolyForMatrix(const BPU_T_GF2_Poly* v) {
	int j, bits_to_print;

	for (j = 0; j < v->elements_in_row; j++) {
		if (j == v->elements_in_row-1) {
			if (v->len % (v->element_bit_size) != 0) {
				bits_to_print = v->len % v->element_bit_size;
			}
			else {
				bits_to_print = v->element_bit_size;
			}
		}
		else {
			bits_to_print = v->element_bit_size;
		}
		BPU_printBinaryLsb(v->elements[j], bits_to_print);
	}
}

void BPU_printGf2Poly(const BPU_T_GF2_Poly* v) {
	int j, bits_to_print;

	fprintf(stderr, "Poly (%4d): ", v->len-1);
	for (j = v->elements_in_row - 1; j >= 0; j--) {
		if (j == v->elements_in_row-1) {
			if (v->len % (v->element_bit_size) != 0) {
				bits_to_print = v->len % v->element_bit_size;
			}
			else {
				bits_to_print = v->element_bit_size;
			}
		}
		else {
			bits_to_print = v->element_bit_size;
		}
		BPU_printBinaryMsb(v->elements[j], bits_to_print);
	}
	fprintf(stderr, "\n");
}

void BPU_printGf2QcMatrix(const BPU_T_GF2_QC_Matrix *v) {
	int ele, i,j,k,s;
	BPU_T_GF2_Poly temp;

	fprintf(stderr, "%s QC Matrix(%i x %i) with %i elements", v->isVertical ? "VERTICAL" : "HORIZONTAL", (v->is_I_appended ? v->k : 0) + v->n, v->k, v->row_element_count * v->column_element_count);
	if (v->is_I_appended)
		fprintf(stderr, " and Identity matrix");
	fprintf(stderr, "\n");

	for (i = 0; i < v->column_element_count; i++) {
		for (k = 0; k < v->element_size; k++) {
			for (j = 0; j < v->row_element_count; j++) {
				int m = j + i * v->row_element_count;
				BPU_gf2PolyCopy(&temp, &v->matrices[m]);
				for (s = 0; s < k; s++) {
					BPU_gf2PolyMulX(&temp);
				}

				BPU_printGf2PolyForMatrix(&temp);
				BPU_gf2PolyFree(&temp, 0);
				fprintf(stderr, "|");
			}

			fprintf(stderr, "\n");
		}

		for (k = 0; k < v->row_element_count * (v->element_size + 1); k++) {
			fprintf(stderr, "-");
		}

		fprintf(stderr, "\n");
	}
}

void BPU_printGf2SparseQcMatrix(const BPU_T_GF2_Sparse_Qc_Matrix *v) {
	int i;
	BPU_T_GF2_Sparse_Poly row;

	fprintf(stderr, "%s QC Matrix(%i x %i) with %i elements\n", v->isVertical ? "VERTICAL" : "HORIZONTAL", v->n, v->k, v->element_count);

	for (i = 0; i < v->k; i++) {
		BPU_gf2SparseQcMatrixGetRow(&row, v, i);
		BPU_printGf2SparsePoly(&row);
		BPU_gf2SparsePolyFree(&row, 0);
	}
}
/* ------------------------------------ Print functions ------------------------------------ */
#endif // BPU_CONF_PRINT

int BPU_gf2VecRand(BPU_T_GF2_Vector *out, int w) {
	int i, j;

	if (w > out->len) {
		BPU_printError("weight error w > l");
		return -2;
	}
	//vector of random weight
	if(w == 0) {
		for(i = 0; i < out->len; i++) {
			BPU_gf2VecSetBit(out, i, BPU_prngGetRand(0, 2));
		}
	}
	//vector of required weight
	else {
		BPU_gf2VecNull(out);

		for(i = 0; i < w; i++) {
			j = BPU_prngGetRand(0, out->len);

			if(BPU_gf2VecGetBit(out, j) == 0) {
				BPU_gf2VecSetBit(out, j, 1);
			}
			else{
				i--;
			}
		}
	}
	return 0;
}

int BPU_gf2VecRandProb(BPU_T_GF2_Vector *out, double prob, int parity) {
	int i, j;
	int hwt = 0;
	double r;

	if (prob > 1 || prob < 0) {
		BPU_printError("probability out of bounds");
		return -2;
	}

	BPU_gf2VecNull(out);
	for (i = 0; i < out->len - 1; i++) {
		r = BPU_prngGetRand(0, RAND_MAX) / (double)RAND_MAX;

		if (r < prob) {
			BPU_gf2VecSetBit(out, i, 1);	
			hwt++;
		}
		else {
			BPU_gf2VecSetBit(out, i, 0);
		}
	}
	
	//setting parity bit
	if (parity == 0) {
		if (hwt % 2 == 0) {
			BPU_gf2VecSetBit(out, out->len - 1, 1);		
			hwt++;
		}
		else
			BPU_gf2VecSetBit(out, out->len - 1, 0);
	}
	else {
		if (hwt % 2 == 0) {
			BPU_gf2VecSetBit(out, out->len - 1, 0);
		}
		else {
			BPU_gf2VecSetBit(out, out->len - 1, 1);
			hwt++;
		}
	}

	return 0;
}

int BPU_gf2MatCopy(BPU_T_GF2_Matrix *out, const BPU_T_GF2_Matrix *in) {
	int i, j;

	if (out->k != in->k || out->n != in->n) {
		BPU_printError("BPU_gf2MatCopy: wrong matrix size");

		return -1;
	}

	// copy the matrix
	for (i = 0; i < in->k; i++) {
		for (j = 0; j < in->elements_in_row; j++) {
			out->elements[i][j] = in->elements[i][j];
		}
	}
	return 0;
}

int BPU_gf2VecPermute(BPU_T_GF2_Vector *vec, const BPU_T_Perm_Vector *permutation) {
	int i;
	BPU_T_GF2_Vector *tmp;

	BPU_gf2VecMalloc(&tmp, vec->len);

	for (i = 0; i < permutation->size; i++) {
		BPU_gf2VecSetBit(tmp, i, BPU_gf2VecGetBit(vec, permutation->elements[i]));
	}
	BPU_gf2VecCopy(vec, tmp);

	BPU_gf2VecFree(&tmp);

	return 0;
}

BPU_T_GF2 BPU_gf2MatGetMaskedBit(const BPU_T_GF2_Matrix *m, uint32_t row, uint32_t bit) {
	int segment, bit_in_segment;

	segment = bit / (m->element_bit_size);
	bit_in_segment = bit % (m->element_bit_size);
	// TODO: consider repleacing di literal 1u
	return m->elements[row][segment] & ((uint32_t)1 << bit_in_segment);
}

BPU_T_GF2 BPU_gf2VecGetMaskedBit(const BPU_T_GF2_Vector *vec, uint32_t bit) {
	int segment, bit_in_segment;

	segment = bit / (vec->element_bit_size);
	bit_in_segment = bit % (vec->element_bit_size);
	// TODO: consider repleacing di literal 1u
	return vec->elements[segment] & ((uint32_t)1 << bit_in_segment);
}

int BPU_gf2MatTransp(BPU_T_GF2_Matrix *out, const BPU_T_GF2_Matrix *in) {
	int i, j;

	if (out->k != in->n || out->n != in->k) {
		BPU_printError("Wrong matrix dimenzion");

		return -1;
	}
	for (i = 0; i < in->k; i++) {
		for (j = 0; j < in->n; j++) {
			BPU_gf2MatSetBit(out, j, i, BPU_gf2MatGetBit(in, i, j));
		} // col loop
	} // row loop

	return 0;
}

void BPU_gf2Swap(BPU_T_GF2 *a, BPU_T_GF2 *b) {
	BPU_T_GF2 tmp;
	tmp = *a;
	*a = *b;
	*b = tmp;
}

void BPU_gf2MatSwapRows(BPU_T_GF2_Matrix *mat, int i, int j) {
	int k;
	for (k = 0; k < mat->elements_in_row; k++) {
		BPU_gf2Swap(&(mat->elements[i][k]), &(mat->elements[j][k]));
	}
}



int BPU_gf2MatFindRow(const BPU_T_GF2_Matrix *mat, int i, int start_index) {
	int k;
	for (k = start_index; k < mat->k; k++) {
		if ( BPU_gf2MatGetMaskedBit(mat, k, i) ){
			return k;
		}
	}
	return -1;
}

int BPU_gf2MatFindCol(const BPU_T_GF2_Matrix *mat, int i, int start_index) {
	int k;
	for (k = start_index; k < mat->n; k++) {
		if ( BPU_gf2MatGetMaskedBit(mat, i, k) ){
			return k;
		}
	}
	return -1;
}

int  BPU_gf2MatMakeSystematic(BPU_T_GF2_Matrix *inout) {
	int act_position = 0;
	int i;
	int row;

	for (act_position = 0; act_position < inout->k; act_position++) {
		row = BPU_gf2MatFindRow(inout, act_position, act_position);
		if (row == -1){
			BPU_printWarning("Not systematic");
			return -1;
		}
		BPU_gf2MatSwapRows(inout, act_position, row);

		// xor with the rest of rows if needed
		for (i = 0; i < inout->k; i++) {
			if ( BPU_gf2MatGetMaskedBit(inout, i, act_position) && act_position != i) {
				BPU_gf2MatXorRows(inout, i, act_position);
			}
		}
	}
	return 0;
}

int BPU_gf2MulMat(BPU_T_GF2_Matrix **out, BPU_T_GF2_Matrix *left, BPU_T_GF2_Matrix *right) {
	int i, j, k;

	if (left->n != right->k) {
		BPU_printError("Incompatible matrix dimensions");
		return -1;
	}

	BPU_gf2MatMalloc(out, left->k, right->n);

	for (i = 0; i < (*out)->k; i++) {
		for (j = 0; j < left->n; j++) {
			if (BPU_gf2MatGetBit(left, i, j)) {
				for (k = 0; k < (*out)->elements_in_row; k++) {
					(*out)->elements[i][k] ^= right->elements[j][k];
				}
			}
		}
	}

	return 0;
}

int BPU_gf2VecConcat(BPU_T_GF2_Vector *out, const BPU_T_GF2_Vector *vec1, const BPU_T_GF2_Vector *vec2) {
	int len = vec1->len + vec2->len;
	int i;

	if (out->len != len) {
		if (BPU_gf2VecResize(out, len)) {
			BPU_printError("resize error");
			return -1;
		}
	}
	else {
		BPU_gf2VecNull(out);
	}
	// copy first vector
	BPU_gf2VecCopy(out, vec1);

	// copy second vector
	for (i = 0; i < vec2->len; i++) {
		BPU_gf2VecSetBit(out, i + vec1->len, BPU_gf2VecGetBit(vec2, i));
	}
	out->len = len;

	return 0;
}

int BPU_gf2VecCrop(BPU_T_GF2_Vector *out, const BPU_T_GF2_Vector *in, const int start, const int length) {
	int i;
	int counter = 0;

	//input test
	if (start < 0 || (length+start) > in->len) {
		BPU_printError("cropVector: bad params");
		return -1;
	}
	if (out->len < length) {
		BPU_printError("cropVector: out vector is smaller then needed");
		return -1;
	}
	for (i = start; i < start+length; i++) {
		BPU_gf2VecSetBit(out, counter, BPU_gf2VecGetBit(in, i));
		counter++;
	}
	return 0;
}

int BPU_gf2MatGetRowAsGf2Vec(BPU_T_GF2_Vector *out, const BPU_T_GF2_Matrix *in, int row) {
	if (out->len != in->n) {
		BPU_printError("dimension is wrong out->len %d != in->n %d", out->len, in->n);
		return -1;
	}
	BPU_gf2MatCopyRowToVec(out, in, row);

	return 0;
}

void BPU_gf2VecCopy(BPU_T_GF2_Vector *dest, const BPU_T_GF2_Vector *src) {
	// if there is not enough space resize
	if (dest->elements_in_row < src->elements_in_row) {
		BPU_gf2VecResize(dest, src->elements_in_row * src->element_bit_size * sizeof(BPU_T_GF2));
	}
	else {
		BPU_gf2VecNull(dest);
	}
	memcpy((void *) (dest->elements), (void *) (src->elements), sizeof(BPU_T_GF2) * (src->elements_in_row));

	dest->len = src->len;
}

int BPU_gf2VecCmp(const BPU_T_GF2_Vector *v1, const BPU_T_GF2_Vector *v2) {
	int i;

	if (v1->len != v2->len) {
		return -1;
	}
	for (i = 0; i < v1->len; i++) {
		if (BPU_gf2VecGetBit(v1, i) != BPU_gf2VecGetBit(v2, i)) {
			return i + 1;
		}
	}
	return 0;
}

int BPU_gf2VecXor(BPU_T_GF2_Vector *out, const BPU_T_GF2_Vector *in) {
	int i;

	if (out->len != in->len) {
		BPU_printError("BPU_gf2VecXor: length error (el. in row) %d != %d, len %d != %d", out->elements_in_row, in->elements_in_row, out->len, in->len);

		return -1;
	}
	for (i = 0; i < out->elements_in_row; i++) {
		out->elements[i] ^= in->elements[i];
	}
	return 0;
}

int BPU_gf2VecMulMat(BPU_T_GF2_Vector *out, const BPU_T_GF2_Vector *v, const BPU_T_GF2_Matrix *b) {
	int i, j;

	if ((v->len != b->k) || (out->len != b->n)) {
		BPU_printError("wrong vector and matrix dimension v->len = %d, b->k = %d", v->len, b->k);
		BPU_printError("wrong vector and matrix dimension out->len = %d, b->n = %d", out->len, b->n);

		return -1;
	}
	// null elements
	BPU_gf2VecNull(out);

	for (i = 0; i < v->len; i++) {
		if (BPU_gf2VecGetBit(v, i)) {
			// xor rows
			for (j = 0; j < out->elements_in_row; j++) {
				out->elements[j] ^= b->elements[i][j];
			}
		}
	}
	return 0;
}

void BPU_gf2MatXorRows(BPU_T_GF2_Matrix *mat, int i, int j) {
	int k;

	for (k = 0; k < mat->elements_in_row; k++) {
		mat->elements[i][k] ^= mat->elements[j][k];
	}
}

int BPU_gf2MatPermute(BPU_T_GF2_Matrix *inout, BPU_T_Perm_Vector *permutation) {
	int i, j;
	int bit = 0;
	BPU_T_GF2_Vector *vector;
	int length;

	if (inout->n != permutation->size) {
		BPU_printError("permutation size not correct m->n = %d, p->size = %d", inout->n, permutation->size);

		return -1;
	}
	length = inout->elements_in_row * (inout->element_bit_size / 8);

	BPU_gf2VecMalloc(&vector, permutation->size);

	for (i = 0; i < inout->k; i++) {
		memcpy((vector->elements), inout->elements[i], length);
		for (j = 0; j < permutation->size; j++) {
			bit = BPU_gf2VecGetBit(vector, permutation->elements[j]);
			BPU_gf2MatSetBit(inout, i, j, bit);
		}
	}
	BPU_gf2VecFree(&vector);
	return 0;
}

int BPU_gf2MatCrop(BPU_T_GF2_Matrix *m, uint16_t width) {
	BPU_T_GF2_Vector *row, *cropped_row;
	int length, i, new_length;

	length = m->elements_in_row * (m->element_bit_size / 8);
	BPU_gf2VecMalloc(&row, m->n);
	BPU_gf2VecMalloc(&cropped_row, width);
	new_length = cropped_row->elements_in_row * (m->element_bit_size / 8);
	for (i = 0; i < m->k; i++) {
		memcpy(row->elements, m->elements[i], length);
		BPU_gf2VecCrop(cropped_row, row, m->k, width);
		free(m->elements[i]);
		m->elements[i] = (BPU_T_GF2*) malloc(new_length);
		memcpy(m->elements[i], cropped_row->elements, new_length);
	}
	m->n = cropped_row->len;
	m->elements_in_row = cropped_row->elements_in_row;

	BPU_gf2VecFree(&row);
	BPU_gf2VecFree(&cropped_row);

	return 0;
}

uint8_t BPU_getParity(BPU_T_GF2 dword) {
	uint32_t tmp = dword;
	tmp = (tmp >> 16) ^ tmp;
	tmp = (tmp >> 8 ) ^ tmp;
	tmp = (tmp >> 4 ) ^ tmp;
	tmp = (tmp >> 2 ) ^ tmp;
	tmp = (tmp >> 1 ) ^ tmp;
	return tmp & 1;
}

void BPU_gf2PolyCopy(BPU_T_GF2_Poly *out, const BPU_T_GF2_Poly *in) {

	int i;
	// allocate output poly
	BPU_gf2PolyMalloc(out, in->len);

	// copy all elements
	if (in->len != 0)
		for (i = 0; i < in->elements_in_row; i++) {
			out->elements[i] = in->elements[i];
		}
}

int BPU_gf2QcMatrixToSparse(BPU_T_GF2_Sparse_Qc_Matrix *out, const BPU_T_GF2_QC_Matrix *in, const int wi[]) {
	int i, counter, bit;

	// allocate output matrix
	BPU_gf2SparseQcMatrixMalloc(out, in->element_count, in->element_size, 1);

	// transpose polynoms
	for (i = 0; i < in->element_count; i++) {
		counter = 0;
		// allocate sparse poly
		BPU_gf2SparsePolyMalloc(&out->matrices[i], wi[i]);
		// set bits
		for (bit = 0; bit < in->matrices[i].len; bit++) {
			if (BPU_gf2VecGetBit(&in->matrices[i], bit) == 1ul) {
				out->matrices[i].index[counter] = (uint32_t)(bit);
				counter++;
			}
		}
		// weight error
		if (counter != wi[i]) {
			BPU_printError("weight error. Weight should be %i, but is %i.\n", wi[i], counter);
			return -1;
		}
	}
	return 0;
}

int BPU_gf2PolyInitRand(BPU_T_GF2_Poly *out, int l, int w, int set_deg) {
	int ret;

	// allocate output poly
	ret = BPU_gf2PolyMalloc(out, l);

	// set random bits
	if (w >= 0)
		ret = BPU_gf2VecRand((BPU_T_GF2_Vector*) out, w);

	// set poly deg
	if (set_deg) BPU_gf2PolySetDeg(out, -1);

	return ret;
}

int BPU_gf2PolyInitRandProb(BPU_T_GF2_Poly *out, int l, double prob, int set_deg, int parity) {
	int ret;

	// allocate output poly
	ret = BPU_gf2PolyMalloc(out, l);

	// set random bits
	if (prob >= 0)
		ret = BPU_gf2VecRandProb((BPU_T_GF2_Vector*)out, prob, parity);

	// set poly deg
	if (set_deg) BPU_gf2PolySetDeg(out, -1);

	return ret;
}

void BPU_gf2SparsePolyCopy(BPU_T_GF2_Sparse_Poly *out, const BPU_T_GF2_Sparse_Poly *in) {
	int i;
	// allocate output poly
	BPU_gf2SparsePolyMalloc(out, in->weight);
	// copy all indexes
	for (i = 0; i < in->weight; i++) {
		out->index[i] = in->index[i];
	}
}

void BPU_gf2SparseQcMatrixTransp(BPU_T_GF2_Sparse_Qc_Matrix *out, const BPU_T_GF2_Sparse_Qc_Matrix *in) {
	int counter, coeff, ele, zero_coeff_is_set;

	// allocate memory for output matrix
	BPU_gf2SparseQcMatrixMalloc(out, in->element_count, in->element_size, 1);

	// for all polynoms
	for (ele = 0; ele < in->element_count; ele++) {
		counter = 0; zero_coeff_is_set = 0;

		// alloc new matrix of same weight
		BPU_gf2SparsePolyMalloc(&out->matrices[ele], in->matrices[ele].weight);

		// is zero coeff set?
		zero_coeff_is_set = (in->matrices[ele].index[0] == 0ul);
		// set zero coeff
		if (zero_coeff_is_set) {
			out->matrices[ele].index[counter] = 0ul;
			counter++;
		}

		// for other coeffs
		for (coeff = in->matrices[ele].weight-1; coeff >= (1-!zero_coeff_is_set); coeff--) {
			out->matrices[ele].index[counter] = in->element_size - in->matrices[ele].index[coeff];
			counter++;
		}
	}
}

void BPU_gf2PolyShiftLeft(BPU_T_GF2_Poly *a, int shift_count) {
	int i;
	int diff, bit_shift, start, shift_right, shift_left;
	uint32_t ele1, ele2;

	// allocate new poly, with additional bits
	// BPU_gf2PolyMalloc(&poly, a->len+shift_count);
	BPU_gf2PolySetDeg(a, a->len+shift_count);

	// calc element shift and bit shift
	diff = shift_count / a->element_bit_size;
	bit_shift = shift_count % a->element_bit_size;

	start = diff;
	// shift elements
	for (i = a->elements_in_row-1; i >= start; i--) {
		// not the first element, concat two elements
		if (i-start != 0) {

			// get first element
			if ((a->element_bit_size - bit_shift) >= a->element_bit_size) {
				ele1 = 0ul;
				shift_right = 0;
			}
			else {
				ele1 = a->elements[i-start-1];
				shift_right = a->element_bit_size - bit_shift;
			}

			// get second element
			if ((i-start) >= a->elements_in_row) {
				ele2 = 0ul;
				shift_left = 0;
			}
			else {
				ele2 = a->elements[i-start];
				shift_left = bit_shift;
			}

			// set element
			a->elements[i] = (ele1 >> shift_right) ^ (ele2 << shift_left);
		}
		// first element, just shift
		else
			a->elements[i] = a->elements[i-start] << bit_shift;
	}

	// set zeros in the beginning
	for (i = 0; i < diff; i++)
		a->elements[i] = 0ul;
}

int BPU_gf2PolyGetHighestBitPos(BPU_T_GF2_Poly *a) {
	int ele;

	// scan all elements and found highest non zero element
	for (ele = a->elements_in_row-1; ele >= 0; ele--) {
		if (a->elements[ele] != 0ul)
			// find highest bit in highest non zero element
			return msb32(a->elements[ele], 1, a->element_bit_size, a->element_bit_size) + ele*a->element_bit_size;
	}

	// poly is zero
	return -1;
}

void BPU_gf2PolySetDeg(BPU_T_GF2_Poly *a, int deg) {
	int j, orig_elements_in_row = a->elements_in_row;

	// find max degree
	if (deg == -1) {
		deg = BPU_gf2PolyGetHighestBitPos(a);
	}

	if (deg != -1) {
		// set degree and element count
		a->len = deg;
		a->elements_in_row = a->len / a->element_bit_size + ((a->len % a->element_bit_size) != 0 ? 1 : 0);
		// reallocate elements
		a->elements = (BPU_T_GF2*) realloc(a->elements, sizeof(BPU_T_GF2) * a->elements_in_row);
		// null new elements
		for (j = orig_elements_in_row; j < a->elements_in_row; j++)
			a->elements[j] = 0ul;
	}
	// poly is zero
	else {
		a->len = 0;
		a->elements_in_row = 0;
	}
}

void BPU_gf2PolyMulX(BPU_T_GF2_Poly *a) {
	int ele;
	uint8_t shift = a->element_bit_size-1;

	// save highest bit
	uint32_t msb = BPU_gf2VecGetBit(a, a->len-1);
	// null highest bit
	BPU_gf2VecSetBit(a, a->len-1, 0ul);

	// for all elements
	for (ele = a->elements_in_row-1; ele >= 1; ele--) {
		a->elements[ele] = (a->elements[ele] << 1) ^ (a->elements[ele-1] >> shift);
	}

	// last element just shift
	a->elements[0] <<= 1;
	// and set lowest bit
	a->elements[0] ^= msb;
}

void BPU_gf2PolyShiftRightOne(BPU_T_GF2_Poly *a) {
	int i;

	// for all elements
	for (i = 0; i < a->elements_in_row-1; i++) {
		// shift right by one and add lowest bit from next element
		a->elements[i] = (a->elements[i] >> 1) ^ ((a->elements[i+1] & 1ul) << (a->element_bit_size-1));
	}

	// last element just shift
	a->elements[a->elements_in_row-1] >>= 1;

}

void BPU_gf2PolyAdd(BPU_T_GF2_Poly *out, const BPU_T_GF2_Poly *in, int crop) {
	int i;

	// if out poly is zero, just copy in into out
	if (out->len == 0) {
		BPU_gf2PolyFree(out, 0);
		BPU_gf2PolyCopy(out, in);
	}
	// if in is non zero
	else if (in->len != 0) {
		// if deg(in) > deg(out)
		if (in->len > out->len)
			// set degree
			BPU_gf2PolySetDeg(out, in->len);
		// make add
		for (i = 0; i < in->elements_in_row; i++)
			out->elements[i] ^= in->elements[i];
	}

	// if set new degree
	if (crop) BPU_gf2PolySetDeg(out, -1);
}

// operation add with binary polynomials with high degree
void BPU_gf2SparsePolyAdd(BPU_T_GF2_Poly *out, const BPU_T_GF2_Sparse_Poly *in) {
	int i;

	// for all coefficients
	for (i = 0; i < in->weight; i++)
		// make add
		BPU_gf2VecSetBit(out, in->index[i], BPU_gf2VecGetBit(out, in->index[i]) ^ 1ul);
}

int BPU_gf2SparsePolyAndHW(const BPU_T_GF2_Poly *a, const BPU_T_GF2_Sparse_Poly *b) {
	int i, hw = 0;

	// for all coefficients
	for (i = 0; i < b->weight; i++)
		// if both poly has set coeff
		if (BPU_gf2VecGetBit(a, b->index[i]) == 1ul)
			// increase hamming weight
			hw++;

	return hw;
}

void BPU_gf2SparsePolyGetShift(BPU_T_GF2_Sparse_Poly *p, const BPU_T_GF2_Sparse_Poly *m, int row_num, int poly_length) {
	int i;

	BPU_gf2SparsePolyCopy(p, m);
	for (i = 0; i < p->weight; i++)
		p->index[i] = ((p->index[i]) + row_num) % poly_length;

}

void BPU_gf2SparseQcMatrixGetRow(BPU_T_GF2_Sparse_Poly *p, const BPU_T_GF2_Sparse_Qc_Matrix *m, int row_num) {
	int i;

	// if rownum is in matrix
	if (row_num < m->k) {
		// allocate output polynom
		BPU_gf2SparsePolyCopy(p, &m->matrices[row_num / m->element_size]);
		// VERTICAL QC matrix
		if (m->isVertical)
			// for all coefficients
			for (i = 0; i < p->weight; i++)
				// shift coefficients
				p->index[i] = ((p->index[i]) + row_num) % m->element_size;
		// horizontal matrix is not supported
		else {
			BPU_printError("BPU_QcMatrixGetRow: HORIZONTAL matrix not supported\n");
		}
	}
	// row num is out of the matrix
	else {
		BPU_printError("BPU_QcMatrixGetRow: row with index %i does not exist\n", row_num);
	}
}

void BPU_gf2PolyMulMod(const BPU_T_GF2_Poly *a, const BPU_T_GF2_Poly *b, BPU_T_GF2_Poly *c, const BPU_T_GF2_Poly *m, int crop) {
	int i;
	BPU_T_GF2_Poly temp_b;

	// if one of factors is zero, product will also be zero
	if (a->len == 0 || b->len == 0) {
		BPU_gf2PolyMalloc(c, 0);
		return;
	}

	// copy b to temp_b and prolong to modulo lenght - 1
	BPU_gf2PolyCopy(&temp_b, b);
	BPU_gf2PolySetDeg(&temp_b, m->len-1);

	// malloc multiplication product
	BPU_gf2PolyMalloc(c, m->len);

	// for length of a
	for (i = 0; i < a->len; i++) {
		// 1 in a -> multiply
		if (BPU_gf2VecGetBit(a, i) == 1ul)
			BPU_gf2PolyAdd(c, &temp_b, 0);

		// multiply b by 2
		BPU_gf2PolyMulX(&temp_b);
	}

	// if do not crop the result length
	if (!crop)
		BPU_gf2PolySetDeg(c, m->len-1);
	else
		BPU_gf2PolySetDeg(c, -1);

	BPU_gf2PolyFree(&temp_b, 0);
}

void BPU_gf2PolyDiv(BPU_T_GF2_Poly *q, BPU_T_GF2_Poly *r, const BPU_T_GF2_Poly *a, const BPU_T_GF2_Poly *b) {
	BPU_T_GF2_Poly divisor, dividend;
	int limit_deg = a->len - b->len;
	int i = 0;

	// copy a, b
	BPU_gf2PolyCopy(&divisor, b);
	BPU_gf2PolyCopy(&dividend, a);

	// allocate quotient
	BPU_gf2PolyMalloc(q, a->len);

	// divisor shift to deg of dividend
	BPU_gf2PolyShiftLeft(&divisor, a->len - b->len);

	while (dividend.len >= b->len && divisor.len > 0) {
		if (dividend.len == divisor.len) {
			BPU_gf2VecSetBit(q, limit_deg - i, 1ul);

			// dividend - divisor
			BPU_gf2PolyAdd(&dividend, &divisor, 1);
		}

		// divisor degree decreased by 1
		BPU_gf2PolyShiftRightOne(&divisor);
		// set actual degree
		BPU_gf2PolySetDeg(&divisor, -1);
		i++;
	}

	BPU_gf2PolyCopy(r, &dividend);
	BPU_gf2PolySetDeg(q, -1);

	// free
	BPU_gf2PolyFree(&dividend, 0);
	BPU_gf2PolyFree(&divisor, 0);
}

// XGCD with binary polynomial with high degree
void BPU_gf2PolyExtEuclidA(BPU_T_GF2_Poly *d, BPU_T_GF2_Poly *s, BPU_T_GF2_Poly *t, const BPU_T_GF2_Poly *a, const BPU_T_GF2_Poly *b, const BPU_T_GF2_Poly *m) {
	BPU_T_GF2_Poly tmp, tmp_2, old_s, old_t, old_r, r, q;
	int deg = (a->len > b->len) ? a->len : b->len;

	// allocate Bezout coeffitients
	BPU_gf2PolyMalloc(s, 0);
	BPU_gf2PolyMalloc(t, deg);
	BPU_gf2PolyMalloc(&old_s, deg);
	BPU_gf2PolyMalloc(&old_t, 0);

	// set initial values
	BPU_gf2PolyCopy(&r, a);
	BPU_gf2PolyCopy(&old_r, b);

	if (a->len == 0) {
		old_t.elements[0] = 1ul;
		BPU_gf2PolySetDeg(&old_t, -1);
	}
	else if (b->len == 0) {
		BPU_gf2PolyFree(&old_r, 0);
		BPU_gf2PolyCopy(&old_r, a);
		old_s.elements[0] = 1ul;
		BPU_gf2PolySetDeg(&old_s, -1);
	}
	// run algoritm, if everything is OK
	else {
		// set initial values
		old_s.elements[0] = 1ul;
		BPU_gf2PolySetDeg(&old_s, -1);
		t->elements[0] = 1ul;
		BPU_gf2PolySetDeg(t, -1);

		BPU_gf2PolySetDeg(&r, -1);
		BPU_gf2PolySetDeg(&old_r, -1);
		// while loop until r is not zero
		while (r.len > 0) {
			// divide
			BPU_gf2PolyDiv(&q, &tmp, &old_r, &r);

			// save old reminder
			BPU_gf2PolyFree(&old_r, 0);
			BPU_gf2PolyCopy(&old_r, &r);

			// save current reminder
			BPU_gf2PolyFree(&r, 0);
			BPU_gf2PolyCopy(&r, &tmp);

			// free
			BPU_gf2PolyFree(&tmp, 0);

			// save s quocient
			BPU_gf2PolyCopy(&tmp, &old_s);
			BPU_gf2PolyFree(&old_s, 0);
			BPU_gf2PolyCopy(&old_s, s);

			BPU_gf2PolyMulMod(&q, s, &tmp_2, m, 1);
			BPU_gf2PolyAdd(&tmp, &tmp_2, 1);
			BPU_gf2PolyFree(s, 0);
			BPU_gf2PolyCopy(s, &tmp);

			// free
			BPU_gf2PolyFree(&tmp, 0);
			BPU_gf2PolyFree(&tmp_2, 0);

			// save t quocient
			BPU_gf2PolyCopy(&tmp, &old_t);
			BPU_gf2PolyFree(&old_t, 0);
			BPU_gf2PolyCopy(&old_t, t);
			BPU_gf2PolyMulMod(&q, t, &tmp_2, m, 1);
			BPU_gf2PolyAdd(&tmp, &tmp_2, 1);
			BPU_gf2PolyFree(t, 0);
			BPU_gf2PolyCopy(t, &tmp);

			// free
			BPU_gf2PolyFree(&tmp_2, 0);
			BPU_gf2PolyFree(&tmp, 0);
			BPU_gf2PolyFree(&q, 0);
		}
	}
	// prepare return values
	BPU_gf2PolyFree(t, 0);
	BPU_gf2PolyFree(s, 0);
	BPU_gf2PolyCopy(d, &old_r);
	BPU_gf2PolyCopy(s, &old_s);
	BPU_gf2PolyCopy(t, &old_t);

	// free
	BPU_gf2PolyFree(&old_s, 0);
	BPU_gf2PolyFree(&old_t, 0);
	BPU_gf2PolyFree(&old_r, 0);
	BPU_gf2PolyFree(&r, 0);
}

int BPU_gf2PolyInv(BPU_T_GF2_Poly *out, const BPU_T_GF2_Poly *a, const BPU_T_GF2_Poly *m) {
	BPU_T_GF2_Poly d, s;
	int ret = 1;

	// call XGCD
	BPU_gf2PolyExtEuclidA(&d, &s, out, a, m, m);

	// if GCD (a,m) is not 1
	if (d.len != 1 || d.elements[0] != 1ul) {
		BPU_printDebug("inverse polynomial NOT found");
		ret = 0;
	}

	// free
	BPU_gf2PolyFree(&d, 0);
	BPU_gf2PolyFree(&s, 0);
	return ret;
}

void BPU_gf2PolyTransp(BPU_T_GF2_Poly *out, const BPU_T_GF2_Poly *in) {
	int i;

	// allocate output poly
	BPU_gf2PolyMalloc(out, in->len);

	// copy zero coefficient
	BPU_gf2VecSetBit(out, 0, BPU_gf2VecGetBit(in, 0));

	// swap other coefficients
	for (i = 1; i < out->len; i++) {
		BPU_gf2VecSetBit(out, out->len-i, BPU_gf2VecGetBit(in, i));
	}
}

void BPU_gf2QcMatrixTransp(BPU_T_GF2_QC_Matrix *out, const BPU_T_GF2_QC_Matrix *in) {
	int i;

	BPU_gf2QcMultirowMatrixMalloc(out, in->column_element_count, in->row_element_count, in->element_size);

	// allocate output matrix
	BPU_gf2QcMatrixMalloc(out, in->row_element_count, in->element_size, !in->isVertical, 0);

	// transpose all polynoms
	for (i = 0; i < in->row_element_count; i++)
		BPU_gf2PolyTransp(&out->matrices[i], &in->matrices[i]);
}

void BPU_gf2QcMultirowMatrixTransp(BPU_T_GF2_QC_Matrix *out, const BPU_T_GF2_QC_Matrix *in) {
	int i, block, position, index;
	int elements_count = in->row_element_count * in->column_element_count;

	BPU_gf2QcMultirowMatrixMalloc(out, in->column_element_count, in->row_element_count, in->element_size);

	for (i = 0; i < elements_count; i++) {
		block = i % in->row_element_count;
		position = i / in->row_element_count;
		index = position + block * in->column_element_count;

		BPU_gf2PolyTransp(&out->matrices[index], &in->matrices[i]);
	}
}

void BPU_gf2QcEyeMatrix(BPU_T_GF2_QC_Matrix *out, int block_size, int blocks_count) {
	int i = 0;
	BPU_gf2QcMultirowMatrixMalloc(out, blocks_count, blocks_count, block_size);

	for (i = 0; i < blocks_count * blocks_count; i++) {
		BPU_gf2PolyMalloc(&out->matrices[i], block_size);
		if (i % (blocks_count + 1) == 0) {
			BPU_gf2VecSetBit(&out->matrices[i], 0, 1);
		}
	}
}

void BPU_gf2QcAppendMatrix(BPU_T_GF2_QC_Matrix *out, BPU_T_GF2_QC_Matrix *first, BPU_T_GF2_QC_Matrix *second, int direction) {
	int i, columns, rows, index, position, block;
	BPU_T_GF2_Poly *matrix;
	if (first->element_size != second->element_size) {
		BPU_printError("BPU_gfQcAppendMatrix: Element size of both matrices must be the same\n");
		return;
	}
	if (direction == 0) { // horizontal
		if (first->column_element_count != second->column_element_count) {
			BPU_printError("BPU_gfQcAppendMatrix: Incompatible dimensions of matrices");
			return;
		}

		rows = first->row_element_count + second->row_element_count;
		columns = first->column_element_count;
	} else { // vertical
		if (first->row_element_count != second->row_element_count) {
			BPU_printError("BPU_gfQcAppendMatrix: Incompatible dimensions of matrices");
			return;
		}

		rows = first->row_element_count;
		columns = first->column_element_count + second->column_element_count;
	}

	BPU_gf2QcMultirowMatrixMalloc(out, rows, columns, first->element_size);

	for (i = 0; i < rows * columns; i++) {
		if (direction == 0) {
			position = i % rows;
			block = i / rows;

			if (position < first->row_element_count) {
				index = position + block * first->row_element_count;
				matrix = &first->matrices[index];
			} else {
				index = (position - first->row_element_count) + block * second->row_element_count;
				matrix = &second->matrices[index];
			}
		} else {
			if (i < first->row_element_count * first->column_element_count) {
				matrix = &first->matrices[i];
			} else {
				index = i - (first->row_element_count * first->column_element_count);
				matrix = &second->matrices[index];
			}
		}

		BPU_gf2PolyCopy(&out->matrices[i], matrix);
	}
}

void BPU_gf2QcSelectMatrix(BPU_T_GF2_QC_Matrix *out, BPU_T_GF2_QC_Matrix *in, int column_start, int column_end, int row_start, int row_end) {
	int i, start_index, row_index, column_index;
	BPU_T_GF2_Poly *temp;

	if (column_start >= column_end || row_start >= row_end) {
		BPU_printError("BPU_gf2QcSelectMatrix: Invalid selection interval\n");
		return;
	}

	if (column_start < 0 || column_end > in->row_element_count || row_start < 0 || row_end > in->column_element_count) {
		BPU_printError("BPU_gf2QcSelectMatrix: Index out of bounds\n");
		return;
	}

	BPU_gf2QcMultirowMatrixMalloc(out, column_end - column_start, row_end - row_start, in->element_size);
	start_index = row_start * in->row_element_count + column_start;

	for (i = 0; i < out->column_element_count * out->row_element_count; i++) {
		row_index = i / out->row_element_count;
		column_index = i % out->row_element_count;
		temp = &in->matrices[start_index + row_index * in->row_element_count + column_index];
		BPU_gf2PolyCopy(&out->matrices[i], temp);
	}
}

void BPU_gf2QcMatrixMultiply(BPU_T_GF2_QC_Matrix *out, BPU_T_GF2_QC_Matrix *left, BPU_T_GF2_QC_Matrix *right) {
	int i, j, row, column;
	BPU_T_GF2_Poly temp, temp_product, module;

	if (left->element_size != right->element_size) {
		BPU_printError("BPU_gf2QcMatrixMultiply: Element size of both matrices must be the same\n");
		return;
	}

	if (left->row_element_count != right->column_element_count) {
		BPU_printError("BPU_gf2QcMatrixMultiply: Incompatible dimensions of matrices");
		return;
	}

	BPU_gf2QcMultirowMatrixMalloc(out, right->row_element_count, left->column_element_count, left->element_size);
	BPU_gf2PolyMalloc(&module, left->element_size + 1);
	BPU_gf2VecSetBit(&module, 0, 1);
	BPU_gf2VecSetBit(&module, left->element_size, 1);

	for (i = 0; i < out->column_element_count * out->row_element_count; i++) {
		BPU_gf2PolyMalloc(&temp, left->element_size);
		row = i / out->row_element_count;
		column = i % out->row_element_count;

		for (j = 0; j < left->row_element_count; j++) {
			BPU_gf2PolyMulMod(&left->matrices[j + row * left->row_element_count], &right->matrices[j * right->row_element_count + column], &temp_product, &module, 0);
			BPU_gf2PolyAdd(&temp, &temp_product, 0);
			BPU_gf2PolyFree(&temp_product, 0);
		}

		BPU_gf2PolyCopy(&out->matrices[i], &temp);
		BPU_gf2PolyFree(&temp, 0);
	}

	BPU_gf2PolyFree(&module, 0);
}

void BPU_gf2QcMatrixAdd(BPU_T_GF2_QC_Matrix *out, BPU_T_GF2_QC_Matrix *in) {
	int i;

	if (out->element_size != in->element_size) {
		BPU_printError("BPU_gf2QcMatrixAdd: Element size of both matrices must be the same\n");
		return;
	}

	if (out->column_element_count != in->column_element_count || out->row_element_count != in->row_element_count) {
		BPU_printError("BPU_gf2QcMatrixAdd: Matrices must have the same dimensions\n");
		return;
	}

	for (i = 0; i < out->column_element_count * out->row_element_count; i++) {
		BPU_gf2PolyAdd(&out->matrices[i], &in->matrices[i], 0);
	}
}

int BPU_gf2QcInvMatrix(BPU_T_GF2_QC_Matrix *out, BPU_T_GF2_QC_Matrix *in) {
	BPU_T_GF2_QC_Matrix eye, appended, rref;
	int found;

	if (in->column_element_count != in->row_element_count) {
		BPU_printError("BPU_gf2QcInvMatrix: Matrix must be a square matrix");
		return 0;
	}

	BPU_gf2QcEyeMatrix(&eye, in->element_size, in->column_element_count);
	BPU_gf2QcAppendMatrix(&appended, in, &eye, 0);
	found = BPU_gf2QcMatrixToRref(&rref, &appended);

	if (found) {
		BPU_gf2QcSelectMatrix(out, &rref, in->column_element_count, in->column_element_count * 2, 0, in->column_element_count);
	}

	BPU_gf2QcMatrixFree(&appended, 0);
	BPU_gf2QcMatrixFree(&rref, 0);
	BPU_gf2QcMatrixFree(&eye, 0);

	return found;
}

int BPU_gf2QcMatrixToRref(BPU_T_GF2_QC_Matrix *out, BPU_T_GF2_QC_Matrix *in) {
	int i, j, k = 0, pivot_found = 0, ret = 1;
	BPU_T_GF2_Poly *temp, pivot, module;
	BPU_gf2QcMatrixCopy(out, in);

	BPU_gf2PolyMalloc(&module, out->element_size + 1);
	BPU_gf2VecSetBit(&module, 0, 1);
	BPU_gf2VecSetBit(&module, out->element_size, 1);

	for (i = 0; i < out->row_element_count; i++) {
		pivot_found = 0;
		for (j = k; j < out->column_element_count; j++) {
			temp = &out->matrices[i + j * out->row_element_count];
			if (!BPU_gf2PolyIsZero(temp)) {
				if (BPU_gf2PolyInv(&pivot, temp, &module) == 1) {
					pivot_found = 1;
					BPU_gf2QcMatrixSwapRowBlocks(out, j, k);
					BPU_gf2QcMatrixMultiplyRow(out, k, &pivot);
					BPU_gf2PolyFree(&pivot, 0);
					break;
				}

				BPU_gf2PolyFree(&pivot, 0);
			}
		}

		if (!pivot_found) {
			if (j == out->column_element_count && k < out->column_element_count) { //this is not satisfying enough
				ret = 0;
				BPU_gf2PolyFree(&module, 0);
				return 0;
			}
			continue;
		}

		for (j = 0; j < in->column_element_count; j++) {
			if (j == k) {
				continue;
			}

			temp = &out->matrices[i + j * out->row_element_count];
			BPU_gf2QcMatrixAddRowMultiple(out, k, j, temp);
		}

		k++;
	}

	BPU_gf2PolyFree(&module, 0);

	return ret;
}

void BPU_gf2QcMatrixSwapRowBlocks(BPU_T_GF2_QC_Matrix *in, int row1, int row2) {
	int i, row1position, row2position;
	BPU_T_GF2_Poly temp;

	if (row1 < 0 || row2 < 0 || row1 > in->column_element_count - 1 || row2 > in->column_element_count - 1) {
		BPU_printError("BPU_gf2QcMatrixSwapRowBlocks: Index out of bounds. \n");
		return;
	}

	if (row1 == row2) {
		return;
	}

	for (i = 0; i < in->row_element_count; i++) {
		row1position = i + row1 * in->row_element_count;
		row2position = i + row2 * in->row_element_count;

		temp = in->matrices[row1position];
		in->matrices[row1position] = in->matrices[row2position];
		in->matrices[row2position] = temp;
	}
}

void BPU_gf2QcMatrixCopy(BPU_T_GF2_QC_Matrix *out, BPU_T_GF2_QC_Matrix *in) {
	int i;
	BPU_gf2QcMultirowMatrixMalloc(out, in->row_element_count, in->column_element_count, in->element_size);
	for (i = 0; i < in->row_element_count * in->column_element_count; i++) {
		BPU_gf2PolyCopy(&out->matrices[i], &in->matrices[i]);
	}
}

void BPU_gf2QcMatrixMultiplyRow(BPU_T_GF2_QC_Matrix *out, int row, BPU_T_GF2_Poly *multiple) {
	int i;
	BPU_T_GF2_Poly temp, module;

	BPU_gf2PolyMalloc(&module, out->element_size + 1);
	BPU_gf2VecSetBit(&module, 0, 1);
	BPU_gf2VecSetBit(&module, out->element_size, 1);

	for (i = 0; i < out->row_element_count; i++) {
		BPU_gf2PolyMulMod(multiple, &out->matrices[i + row * out->row_element_count], &temp, &module, 0);
		BPU_gf2PolyFree(&out->matrices[i + row * out->row_element_count], 0);
		BPU_gf2PolyCopy(&out->matrices[i + row * out->row_element_count], &temp);
		BPU_gf2PolyFree(&temp, 0);
	}
	BPU_gf2PolyFree(&module,0);
}

void BPU_gf2QcMatrixAddRowMultiple(BPU_T_GF2_QC_Matrix *out, int source, int target, BPU_T_GF2_Poly *multiple) {
	int i;
	BPU_T_GF2_Poly temp, tmp_mul, module;

	BPU_gf2PolyCopy(&tmp_mul, multiple);
	BPU_gf2PolyMalloc(&module, out->element_size + 1);
	BPU_gf2VecSetBit(&module, 0, 1);
	BPU_gf2VecSetBit(&module, out->element_size, 1);

	for (i = 0; i < out->row_element_count; i++) {
		BPU_gf2PolyMulMod(&tmp_mul, &out->matrices[i + source * out->row_element_count], &temp, &module, 0);
		BPU_gf2PolyAdd(&out->matrices[i + target * out->row_element_count], &temp, 0);
		BPU_gf2PolyFree(&temp, 0);
	}

	BPU_gf2PolyFree(&module, 0);
	BPU_gf2PolyFree(&tmp_mul, 0);
}

void BPU_gf2QcPermuteMatrixColumns(BPU_T_GF2_QC_Matrix *out, BPU_T_Perm_Vector *permutation) {
	int* permutation_track = (int*)calloc(permutation->size, sizeof(int));
	int i, j, p, t;
	BPU_T_GF2_Poly temp;
	for (i = 0; i < permutation->size; i++) {
		permutation_track[i] = i;
	}

	for (i = 0; i < permutation->size; i++) {
		for (j = 0; j < permutation->size; j++) {
			if (permutation_track[j] == permutation->elements[i]) {
				p = j;
				break;
			}
		}

		t = permutation_track[i];
		permutation_track[i] = permutation_track[p];
		permutation_track[p] = t;

		for (j = 0; j < out->column_element_count; j++) {
			temp = out->matrices[j * out->row_element_count + i];
			out->matrices[j * out->row_element_count + i] = out->matrices[j * out->row_element_count + p];
			out->matrices[j * out->row_element_count + p] = temp;
		}
	}
	free(permutation_track);
}

void BPU_gf2QcPermuteMatrixRows(BPU_T_GF2_QC_Matrix *out, BPU_T_Perm_Vector *permutation) {
	int* permutation_track = (int*)calloc(permutation->size, sizeof(int));
	int i, j, p, t;
	BPU_T_GF2_Poly temp;
	for (i = 0; i < permutation->size; i++) { 
		permutation_track[i] = i;
	}
	for (i = 0; i < permutation->size; i++) {
		for (j = 0; j < permutation->size; j++) {
			if (permutation_track[j] == permutation->elements[i]) {
				p = j;
				break;
			}
		}

		t = permutation_track[i];
		permutation_track[i] = permutation_track[p];
		permutation_track[p] = t;

		for (j = 0; j < out->row_element_count; j++) {
			temp = out->matrices[j + out->column_element_count * i];
			out->matrices[j + out->column_element_count * i] = out->matrices[j + out->column_element_count * p];
			out->matrices[j + out->column_element_count * p] = temp;
		}
	}
	free(permutation_track);
}

int BPU_gf2PolyIsZero(const BPU_T_GF2_Poly *a) {
	int i;

	// scan all elements
	for (i = 0; i < a->elements_in_row; i++)
		// if there is non zero element, poly is not zero
		if (a->elements[i] != 0ul)
			return 0;

	// all elements are zero, so also poly is zero
	return 1;
}
