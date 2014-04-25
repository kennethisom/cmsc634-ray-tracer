The application implements all functionality described in the project description.  To use the depth of field capabilities provide a floating point number using the -a flag and the number of rays per pixel using the -r flag.  Example:

trace -a 0.1 -r 50 < input.nff

There is a slight issue with the coloring of the circle that has transparency in the test.nff file (compared to the provided example).  But transparency seems to be working correctly on gears2.nff.

I used the following sites to brush up on my C++:

http://www.cplusplus.com/doc/tutorial/classes/
http://www.cplusplus.com/reference/vector/vector/
http://www.cplusplus.com/reference/vector/vector/vector/
http://www.cplusplus.com/doc/tutorial/pointers/
http://www.cplusplus.com/reference/algorithm/max/
http://msdn.microsoft.com/en-us/library/17w5ykft.aspx
http://stackoverflow.com/questions/18494218/how-to-convert-char-into-float
http://www.cplusplus.com/reference/cstdlib/exit/
http://stackoverflow.com/questions/686353/c-random-float-number-generation
http://www.cplusplus.com/reference/cstdlib/rand/
http://en.wikipedia.org/wiki/Unit_circle
http://www.cplusplus.com/reference/cmath/sin/