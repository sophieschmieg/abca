package fields;

import java.util.List;

public interface Ring<T extends Element> extends MathSet<T> {
	public T zero();
	public T one();
	public int characteristic();
	public T add(T t1, T t2);
	public T negative(T t);
	public T multiply(T t1, T t2);
	public boolean isUnit(T t);
	public int getNumberOfUnits();
	public T inverse(T t);
	public boolean isIntegral();
	public boolean isEuclidean();
	public List<T> quotientAndRemainder(T dividend, T divisor);
	public Group<T> getAdditiveGroup();
	public Group<T> getMultiplicativeGroup();
	public T add(T t1, T t2, T t3);
	public T subtract(T minuend, T subtrahend);
	public T getInteger(int n);
	public T multiply(int n, T t);
	public T multiply(T t1, T t2, T t3);
	public T multiply(int n, T t1, T t2);
	public T multiply(int n, T t1, T t2, T t3);
	public T power(T t, int n);
	public T gcd(T t1, T t2);
	public List<T> extendedEuclidean(T t1, T t2);
	public Iterable<T> getUnits();
}
