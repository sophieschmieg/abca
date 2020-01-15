package fields;

public interface Group<T extends Element> extends MathSet<T> {
	public T neutral();
	public T inverse(T t);
	public T operate(T t1, T t2);
}
