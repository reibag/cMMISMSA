/**
 *
 *      @file masks.c
 *      @brief Atom selector masks
 *
 *      @author Alvaro Cortes Cabrera <alvarocortesc@gmail.com<
 *      @date 2021/03
 *
 *      This program is free software; you can redistribute it and/or modify
 *      it under the terms of the GNU General Public License as published by
 *      the Free Software Foundation version 2 of the License.
 *
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU General Public License for more details.
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void select_on_mask(MOL2 **mymol, char *mask, int verbose)
{
	MOL2 *mol = NULL;
	char *token = NULL;
	char *mask_copy = NULL;
	char *dash = NULL;
	int first = -1, last = -1, i = 0;
	int total = 0;

	if (mymol == NULL || *mymol == NULL)
		return;

	mol = *mymol;
	if (mask == NULL || *mask == '\0')
	{
		for (i = 0; i < mol->n_atoms; i++)
		{
			mol->backbone[i] = 1;
			++total;
		}
		if (verbose)
		{
			fprintf(stderr, "Mask: Selected %i atoms (default all atoms)\n", total);
			fflush(stderr);
		}
		return;
	}

	mask_copy = strdup(mask);
	if (mask_copy == NULL)
		return;

	token = strtok(mask_copy, ",");
	while (token != NULL)
	{
		first = -1;
		last = -1;
		dash = strchr(token, '-');
		if (dash != NULL)
		{
			*dash = '\0';
			first = atoi(token);
			last = atoi(dash + 1);
		}
		else
		{
			first = atoi(token);
			last = first;
		}

		if (last == -1)
			last = first;

		for (i = 0; i < mol->n_atoms; i++)
		{
			if (mol->internal_res_num[i] + 1 >= first && mol->internal_res_num[i] + 1 <= last)
			{
				mol->backbone[i] = 1;
				++total;
			}
		}

		token = strtok(NULL, ",");
	}

	free(mask_copy);

	if (verbose)
	{
		fprintf(stderr, "Mask: Selected %i atoms\n", total);
		fflush(stderr);
	}
}

void select_on_mask_atoms(MOL2 **mymol, char *mask, int verbose)
{
	MOL2 *mol = NULL;
	char *token = NULL;
	char *mask_copy = NULL;
	char *dash = NULL;
	int first = -1, last = -1, i = 0;
	int total = 0;

	if (mymol == NULL || *mymol == NULL)
		return;

	mol = *mymol;
	if (mask == NULL || *mask == '\0')
		return;

	mask_copy = strdup(mask);
	if (mask_copy == NULL)
		return;

	token = strtok(mask_copy, ",");
	while (token != NULL)
	{
		first = -1;
		last = -1;
		dash = strchr(token, '-');
		if (dash != NULL)
		{
			*dash = '\0';
			first = atoi(token);
			last = atoi(dash + 1);
		}
		else
		{
			first = atoi(token);
			last = first;
		}

		if (last == -1)
			last = first;
		if (last >= mol->n_atoms)
			last = mol->n_atoms - 1;
		if (first < 1)
			first = 1;

		for (i = first - 1; i <= last - 1; i++)
		{
			if (i >= 0 && i < mol->n_atoms)
			{
				mol->backbone[i] = 1;
				++total;
			}
		}

		token = strtok(NULL, ",");
	}

	free(mask_copy);

	if (verbose)
	{
		fprintf(stderr, "Mask: Selected %i atoms\n", total);
		fflush(stderr);
	}
}

