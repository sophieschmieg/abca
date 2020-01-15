package fields;


public interface Field<T extends Element> extends MathSet<T> {
	public T zero();
	public T one();
	public int characteristic();
	public T add(T t1, T t2);
	public T negative(T t);
	public T multiply(T t1, T t2);
	public T inverse(T t);
	public Group<T> getAdditiveGroup();
	public Group<T> getMultiplicativeGroup();
	public T add(T t1, T t2, T t3);
	public T subtract(T minuend, T subtrahend);
	public T getInteger(int n);
	public T multiply(int n, T t);
	public T multiply(T t1, T t2, T t3);
	public T multiply(int n, T t1, T t2);
	public T multiply(int n, T t1, T t2, T t3);
	public T divide(T dividend, T divisor);
	public T power(T t, int n);
	public Iterable<T> getNonZeroElements() throws InfinityException;
	public boolean hasRoot(T t, int n);
	public T root(T t, int n);
        public boolean hasSqrt(T t);
        public T sqrt(T t);
	public T primitiveRoot();
}
