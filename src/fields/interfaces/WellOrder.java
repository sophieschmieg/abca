package fields.interfaces;

public interface WellOrder<T extends Element<T>> extends MathSet<T> {
	public T least();
	public T increment(T t);
}
