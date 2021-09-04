package fields.quaternions;

import java.util.List;

import fields.interfaces.Algebra;
import fields.interfaces.Element;
import fields.interfaces.MathMap;
import fields.interfaces.VectorSpace;
import fields.quaternions.AbstractQuaternions.Quaternion;
import fields.vectors.Matrix;
import fields.vectors.SubVectorSpace;

public interface Quaternions<T extends Element<T>> extends Algebra<T, Quaternion<T>>, VectorSpace<T, Quaternion<T>>{
	Quaternion<T> i();
	Quaternion<T> j();
	Quaternion<T> k();

	T a();
	T b();

	T discriminant();
	int hilbertSymbol();
	boolean isSplit();

	Quaternion<T> getElement(T t, T x, T y, T z);
	T asRealPart(Quaternion<T> t);

	Quaternion<T> conjugate(Quaternion<T> t);
	
	SubVectorSpace<T, Quaternion<T>> pureQuaternions();

	T reducedNorm(Quaternion<T> t);
	T reducedTrace(Quaternion<T> t);
	
	QuadraticForm<T, Quaternion<T>> normForm();
	QuadraticForm<T, Quaternion<T>> restrictedNormForm();

	List<Matrix<T>> matrixBasis();
	Matrix<T> asMatrix(Quaternion<T> t);
	MathMap<Quaternion<T>, Matrix<T>> asMatrixMap();
}
