## Copyright  2009-2020  Marc Hofmann and Achim Zeileis
##
## This file is part of the 'lmSubsets' R extension.
##
## 'lmSubsets' is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## 'lmSubsets' is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with 'lmSubsets'.  If not, see <http://www.gnu.org/licenses/>.



lmSubsets <- function (formula, ...)
    UseMethod("lmSubsets")

lmSelect <- function (formula, ...)
    UseMethod("lmSelect")


refit <- function (object, ...)
    UseMethod("refit")


model_response <- function (data, ...)
    UseMethod("model_response")

model_response.default <- function (data, type = "any", ...)
    stats::model.response(data = data, type = type)
