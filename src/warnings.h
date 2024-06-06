/*
 * Copyright (C) 2019  Leo Singer
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 * Portable preprocessor macros for line-level compiler warning control.
 */

#ifndef WARNINGS_H
#define WARNINGS_H

#ifdef __ICC
#define WARNINGS_PUSH _Pragma("warning(push)")
#define WARNINGS_POP _Pragma("warning(pop)")
#else  /* gcc or clang */
#define WARNINGS_PUSH _Pragma("GCC diagnostic push")
#define WARNINGS_POP _Pragma("GCC diagnostic pop")
#endif

#ifdef __ICC
#define WARNINGS_IGNORE_INCOMPATIBLE_POINTER_TYPES _Pragma("warning(disable:167)")
#else  /* gcc or clang */
#define WARNINGS_IGNORE_INCOMPATIBLE_POINTER_TYPES _Pragma("GCC diagnostic ignored \"-Wincompatible-pointer-types\"")
#endif

#ifdef __ICC
#define WARNINGS_IGNORE_CAST_ALIGN
#else  /* gcc or clang */
#define WARNINGS_IGNORE_CAST_ALIGN _Pragma("GCC diagnostic ignored \"-Wcast-align\"")
#endif

#endif /* WARNINGS_H */
