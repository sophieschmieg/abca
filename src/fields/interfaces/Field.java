package fields.interfaces;

import java.math.BigInteger;
import java.util.List;

import fields.exceptions.InfinityException;

public interface Field<T extends Element<T>> extends MathSet<T>, Ring<T> {
	public T zero();
	public T one();
	public BigInteger characteristic();
	public T add(T t1, T t2);
	public T negative(T t);
	public T multiply(T t1, T t2);
	public T inverse(T t);
	public Group<T> getAdditiveGroup();
	public Group<T> getMultiplicativeGroup();
	public T add(T t1, T t2, T t3);
	public T subtract(T minuend, T subtrahend);
	public T getInteger(int n);
	public T getInteger(BigInteger n);
	public T multiply(int n, T t);
	public T multiply(BigInteger n, T t);
	public T multiply(T t1, T t2, T t3);
	public T multiply(int n, T t1, T t2);
	public T multiply(BigInteger n, T t1, T t2);
	public T multiply(int n, T t1, T t2, T t3);
	public T multiply(BigInteger n, T t1, T t2, T t3);
	public T divide(T dividend, T divisor);
	public T power(T t, int n);
	public T power(T t, BigInteger n);
	public Iterable<T> getNonZeroElements() throws InfinityException;
	public List<Polynomial<T>> factorization(Polynomial<T> t);
	public boolean hasRoots(Polynomial<T> t);
	public List<T> roots(Polynomial<T> t);
	public boolean hasRoot(T t, int n);
	public List<T> roots(T t, int n);
	public boolean hasSqrt(T t);
	public List<T> sqrt(T t);
	public T primitiveRoot();
}
