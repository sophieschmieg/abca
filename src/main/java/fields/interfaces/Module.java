package fields.interfaces;

import java.math.BigInteger;
import java.util.List;

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
	public List<T> nonTrivialCombination(List<S> s);
	public List<List<T>> nonTrivialCombinations(List<S> s);
	public List<S> getModuleGenerators();
	public default List<List<T>> getModuleGeneratorRelations() {
		return nonTrivialCombinations(getModuleGenerators());
	}
	public Vector<T> asVector(S s);
	public S fromVector(Vector<T> vector);
}
