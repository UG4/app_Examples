-- Copyright (c) 2010-2016:  G-CSC, Goethe University Frankfurt
-- Authors: Sebastian Reiter
-- 
-- This file is part of UG4.
-- 
-- UG4 is free software: you can redistribute it and/or modify it under the
-- terms of the GNU Lesser General Public License version 3 (as published by the
-- Free Software Foundation) with the following additional attribution
-- requirements (according to LGPL/GPL v3 §7):
-- 
-- (1) The following notice must be displayed in the Appropriate Legal Notices
-- of covered and combined works: "Based on UG4 (www.ug4.org/license)".
-- 
-- (2) The following notice must be displayed at a prominent place in the
-- terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
-- 
-- (3) The following bibliography is recommended for citation and must be
-- preserved in all covered files:
-- "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
--   parallel geometric multigrid solver on hierarchically distributed grids.
--   Computing and visualization in science 16, 4 (2013), 151-164"
-- "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
--   flexible software system for simulating pde based models on high performance
--   computers. Computing and visualization in science 16, 4 (2013), 165-179"
-- 
-- This program is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
-- GNU Lesser General Public License for more details.


--------------------------------------------------------------------------------
-- Section 0: comments
--------------------------------------------------------------------------------
-- Lines starting with a -- are regarded as comments in lua.
-- A comment does not have to start at the beginning of a line but may be
-- appended to a regular code statement.
-- there are also
--[[ 
multi-line
comments
]]--


--------------------------------------------------------------------------------
print("SECTION 1: Variable types and type conversion")

-- lua variables do not have a fixed type. A variable is declared by simply using it:
var = 3
print(var)
var = "now var is a string"
print(var)

-- One can use the .. operator to concatenate strings
-- Note that known types are implicitly converted to strings on the fly
var1 = "A concatenated"
var2 = " message. And 5 = "
var3 = 5
str = var1 .. var2 .. var3
print(str)

-- a string containing a number can be converted to a number using 'tonumber':
a = "5"
b = tonumber(a)
print(a .. "+" .. a .. " = " .. b + b)



--------------------------------------------------------------------------------
print("")
print("SECTION 2: Control structures")
-- lua supports the common control structures like for, if etc.

-- A check for equality is done with ==
if varUninitialized == nil then
	print("Uninitalized variables default to 'nil'.")
else 
	print("varUninitialized isn't uninitialized!")
end

-- the negation operator in lua is ~
if 5 ~= 4 then
	print("5 does not equal 4")
end

-- We now use a for loop to assemble a string.
-- We iterate from 1 to 5. The default step size is 1.
str = ""
for i = 1, 5 do
	str = str .. i .. " "
end
print("assembled string in for loop: " .. str)

-- Now lets use a step size of 2
str = ""
for i = 1, 5, 2 do
	str = str .. i .. " "
end
print("assembled sting in for loop with step size 2: " .. str)

-- There are many more control structures in lua. Please check the lua reference
-- for more details.



--------------------------------------------------------------------------------
print("")
print("SECTION 3: Tables")

-- A simple array. Note that indices start at 1.
write("An array: ") -- in contrary to 'print', write doesn't automatically do line breaks
a = {-1, -2, -3, -4, -5}
i = 1
while a[i] ~= nil do
	write(a[i], ", ")
	i = i + 1
end
write("\n") -- newline (line break)


-- Arrays can also be filled dynamically. Note that the variable first has to
-- contain a table.
-- Another interesting thing about lua tables (and thus arrays) is that different
-- entries can contain values of different types.
dynArray = {}
dynArray[1] = "a string "
dynArray[2] = "and a number: "
dynArray[3] = 13

-- we're using the ipairs method here, which allows us to iterate over all entries
-- of an array until the value nil is reached.
str = ""
for i, tval in ipairs(dynArray) do
	str = str .. tval
end

print("Concatenated dynArray: " .. str)

-- tables can of course also contain other tables:
dynArray[4] = {33, 35, 1100, "end of subtable"}


-- Note that arrays are not copied during assignment. Instead both variables
-- then point to the same table
tmpArray = dynArray
-- tmpArray and dynArray now both point to the same table.


-- tables can also be used as associative containers (indeed, this is what they are):
con = {}
con["name"] = "Somebody"
con["occupation"] = "Something"
con["age"] = 99

-- There's an interesting alternative syntax to access those fields
-- which are identified by a string:
print("Name: " .. con.name)
print("Occupation: " .. con.occupation)
-- of course the original array-like access is also valid:
print("Age: " .. con["age"])


-- Be sure to check the reference documentation on tables.


--------------------------------------------------------------------------------
print()
print("SECTION 4: Functions")

-- Functions are declared by the 'function' key word:
function AddFunc(a, b)
	return a + b
end

sum = AddFunc(48, 89)
print("48 and 89 sums up to " .. sum)


-- Lua functions can have an arbitrary number of return values. Return values
-- are separated by commas.
function MultTwoVars(x, y, a)
	return x*a, y*a
end

x, y = MultTwoVars(3, 4, 2)
print("3 and 4 multiplied by 2 results in " .. x .. " and " .. y)


-- omitted parameters are set to nil:
function myprint(a, b)
	if(a ~= nil) then print(a) end
	if(b ~= nil) then print(b) end	
end

myprint("omitting parameter 2 in 'myprint'") -- omitted parameter 'b'

-- (this is also true for wrapped C/C++-functions)

-- functions can also be part of tables (this can be used to create 'namespaces')
-- The following syntax makes sure that a table is not overwritten if it was
-- already defined somewhere else. This works, because 'or' returns the first
-- expression which evaluates to true.
namespace = namespace or {} 
function namespace.print_hello()
	print("namespace prints hello")
end
namespace.print_hello();


--------------------------------------------------------------------------------
print()
print("SECTION 5: Scopes")
-- Lua's scope rules differ from many other (classical) programming languages.
-- A variable used in a method is by default a global variable and exists even
-- after the method was terminated.
-- If you want to use temporary variables only visible inside a function
-- itself, you have to use the 'local' keyword (recommended!):

-- An example with a non-local variable first used in a function:
function SetGlobalVar(a)
	globalVar = a
end

print("Value of globalVar before call to SetGlobalVar(77): ")
print(globalVar)
SetGlobalVar(77)
print("Value of globalVar after call to SetGlobalVar(77): ")
print(globalVar)


-- And now an example with a local variable
print()
function SetLocalVar(a)
	local localVar = a
end

print("Value of localVar before call to SetLocalVar(33): ")
print(localVar)
SetLocalVar(33)
print("Value of localVar after call to SetLocalVar(33): ")
print(localVar)


--------------------------------------------------------------------------------
print()
print("SECTION 6: Standard libraries")

-- string manipulation (http://www.lua.org/manual/5.1/manual.html#5.4)
-- sometimes you want to format your strings like this
print("normal: " .. math.pi)
print("formatted: " .. string.format ("%.40f", math.pi) )

-- math functions (http://www.lua.org/manual/5.1/manual.html#5.6)
print("math functions:")
print("math.abs     math.acos    math.asin    math.atan    math.atan2")
print("math.ceil    math.cos     math.deg     math.exp     math.floor")
print("math.log     math.log10   math.max     math.min     math.mod")
print("math.pow     math.rad     math.sin     math.sqrt    math.tan")
print("math.frexp   math.ldexp   math.random  math.randomseed")

-- i/o . (http://www.lua.org/manual/5.1/manual.html#5.7)
-- with io.open, you get access to files. io.open returns a handle which is actually a FILE* pointer
print("writing files ''")
file = io.open("sample-file-from-lua-programming.txt", "a")
-- use file:write to write data to the file
file:write("hello world")  -- this is the same as file:write(tostring(s))
io.close(file)

-- operating system facilities (http://www.lua.org/manual/5.1/manual.html#5.8)
-- we only cover time and date functions here:
print("clock example with string concatenation.")
tBefore = os.clock()
-- do something
for i = 1, 10000, 1 do
str = str .. i .. " "
end
tAfter = os.clock()
print("Repeated concatenation took " .. tAfter-tBefore .. " seconds")

print("it is " .. os.date("%X on %x"))
