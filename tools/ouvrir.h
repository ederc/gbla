/* gbla: Gröbner Basis Linear Algebra
 * Copyright (C) 2015 Brice Boyer <brice.boyer@lip6.fr>
 *                    Jean-Charles Faugère <jean-charles.faugere@inria.fr>
 * This file is part of gbla.
 * gbla is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * gbla is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with gbla . If not, see <http://www.gnu.org/licenses/>.
 */


#ifndef __GBLA_ouvrir_H
#define __GBLA_ouvrir_H

FILE* ouvrir(const char* a,char* b)
{
	if ( (strcmp(a,"-") == 0) && (strcmp(b,"r") == 0) )
		return stdin;
	else
	{
		FILE *f=fopen(a,b);
		if (!f)
		{
			fprintf(stderr,"Can't open %s\n",a);
			exit(1);
		}
		return f;
	}
}

#endif

/* vim: set ft=c: */
