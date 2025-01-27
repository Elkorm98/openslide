/*
 *  OpenSlide, a library for reading whole slide image files
 *
 *  Copyright (c) 2014 Carnegie Mellon University
 *  All rights reserved.
 *
 *  OpenSlide is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, version 2.1.
 *
 *  OpenSlide is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with OpenSlide. If not, see
 *  <http://www.gnu.org/licenses/>.
 *
 */

#ifdef _WIN32
#define _WIN32_WINNT 0x0600

#include <stdio.h>
#include <windows.h>

// avoid warning about missing prototype
int wmain(int argc, wchar_t **argv);

int wmain(int argc, wchar_t **argv) {
  if (argc != 3) {
    fwprintf(stderr, L"Usage: %s <src> <dst>\n", argv[0]);
    return 1;
  }

  const wchar_t *src = argv[1];
  const wchar_t *dst = argv[2];
  if (!CreateSymbolicLink(dst, src, 0)) {
    fprintf(stderr, "Failed with error %lu\n", GetLastError());
    return 1;
  }
  return 0;
}

#else

int main ()
{
  return 0;
}

#endif
