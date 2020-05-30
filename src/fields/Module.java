package fields;

import java.math.BigInteger;

public interface Module<T extends Element, S extends Element> extends MathSet<S> {
	public Ring<T> getRing();
	public S zero();
	public S add(S s1, S s2);
	public S negative(S s);
	public S scalarMultiply(T t, S s);
	public boolean isFree();
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
}
