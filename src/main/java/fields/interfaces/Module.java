package fields.interfaces;

import java.math.BigInteger;
import java.util.List;

import fields.vectors.Matrix;
import fields.vectors.MatrixModule;
import fields.vectors.Vector;

public interface Module<T extends Element<T>, S extends Element<S>> extends MathSet<S> {
	public Ring<T> getRing();
	public S zero();
	public S add(S s1, S s2);
	public S negative(S s);
	public S scalarMultiply(T t, S s);
	public boolean isFree();
	public boolean isTorsionFree();
	public Ideal<T> annihilator();
	public Group<S> getAdditiveGroup();
	public S add(S s1, S s2, S s3);
	public S subtract(S minuend, S subtrahend);
	public S scalarMultiply(int n, S s);
	public S scalarMultiply(BigInteger n, S s);
	public S scalarMultiply(T t1, T t2, S s);
	public S scalarMultiply(int n, T t, S s);
	public S scalarMultiply(BigInteger n, T t, S s);
	public S scalarMultiply(int n, T t1, T t2, S s);
	public S scalarMultiply(BigInteger n, T t1, T t2, S s);
	public boolean isLinearIndependent(List<S> s);
	public boolean isGeneratingModule(List<S> s);
	public List<Vector<T>> nonTrivialCombinations(List<S> s);
	public boolean isSubModuleMemberModule(MatrixModule<T> module, Matrix<T> m, Vector<T> b);
	public Vector<T> asSubModuleMemberModule(MatrixModule<T> module, Matrix<T> m, Vector<T> b);
	public List<Vector<T>> syzygyProblemModule(MatrixModule<T> module, Matrix<T> m);
	public List<Vector<T>> simplifySubModuleGeneratorsModule(MatrixModule<T> module, Matrix<T> m);
	public boolean isSubModuleMemberModule(List<S> m, S b);
	public S asSubModuleMemberModule(List<S> m, S b);
	public List<Vector<T>> syzygyProblemModule(List<S> m);
	public List<S> simplifySubModuleGeneratorsModule(List<S> m);
	public List<S> getModuleGenerators();
	public List<Vector<T>> getSyzygies();
	public Vector<T> asVector(S s);
	public S fromVector(Vector<T> vector);
}
