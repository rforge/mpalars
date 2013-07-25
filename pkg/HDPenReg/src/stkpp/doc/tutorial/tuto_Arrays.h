/* kennel tutorial: how to used the STK++ kernel tools
 **/
namespace STK
{

/** @page PageArrays The STK++ Arrays
 *
 * The containers/arrays you use in order to store and process the data in your
 * application greatly influence the speed and the memory usage of your application.
 * STK++ propose a large choice of containers/arrays and methods that you can used in
 * conjunction with them.
 *
 * @p There are mainly two kinds of arrays you can use with STK++:
 * @li The %Array2D classes which are the classes defined in the oldest versions of STK++,
 * @li the %CArray class which have been introduced in version 0.4 of STK++ library.
 *
 * Before explaining the usage and differences between the different arrays,
 * we first introduce some vocabulary.
 * The terminology used in STK++ project for the arrays are the following:
 * vectors are just a special case of matrices, with either 1 row or 1 column.
 * @li In the case where they have 1 column, such vectors are called column-vectors,
 * often abbreviated just as @em vectors,
 * @li in the other case, where they have 1 row, they are called row-vectors, often
 * abbreviated just as @em points.
 *
 * The %Array2D classes are very flexible if you need to add, insert, remove, resize,...
 * quickly rows or columns to your container. On the other hand the %CArray classes
 * allow you to perform fast and inlined operations (whenever possible) on matrices
 * and vectors. Moreover, the storing scheme allow you to interface them easily
 * to other linear algebra libraries (e.g. Lapack, Blas, ...).
 *
 * <li> @ref PageArrays
 * <ul>
 * <li> @ref ArraysFirstExemple
 * <li>
 * </ul>
 *
 **/

} // namespace STK
