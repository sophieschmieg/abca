package fields.interfaces;

import java.math.BigInteger;
import java.util.List;

import fields.polynomials.UnivariatePolynomialRing;

public interface Ring<T extends Element<T>> extends MathSet<T> {
	public T zero();
	public T one();
	public BigInteger characteristic();
	public T add(T t1, T t2);
	public T negative(T t);
	public T multiply(T t1, T t2);
	public boolean isUnit(T t);
	public BigInteger getNumberOfUnits();
	public T inverse(T t);
	public boolean isIntegral();
	public boolean isZeroDivisor(T t);
	public boolean isEuclidean();
	public boolean isDivisible(T dividend, T divisor);
	public List<T> quotientAndRemainder(T dividend, T divisor);
	public T divide(T dividend, T divisor);
	public T remainder(T dividend, T divisor);
	public BigInteger euclidMeasure(T t);
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
	public T power(T t, int n);
	public T power(T t, BigInteger n);
	public T projectToUnit(T t);
	public T upToUnit(T t);
	public T gcd(T t1, T t2);
	public T lcm(T t1, T t2);
	public List<T> extendedEuclidean(T t1, T t2);
	public Iterable<T> getUnits();
	public int krullDimension();
	public Ideal<T> getIdeal(List<T> generators);
	public Ideal<T> add(Ideal<T> t1, Ideal<T> t2);
	public Ideal<T> multiply(Ideal<T> t1, Ideal<T> t2);
	public Ideal<T> intersect(Ideal<T> t1, Ideal<T> t2);
	public UnivariatePolynomialRing<T> getUnivariatePolynomialRing();
}
