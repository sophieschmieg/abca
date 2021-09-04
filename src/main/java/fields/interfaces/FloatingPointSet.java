package fields.interfaces;

public interface FloatingPointSet<T extends Element<T>, S extends MathSet<T>> extends MathSet<T> {
	boolean close(T t1, T t2);
	int precision();
	T roundToFixedPoint(T t, int precision);
	S withPrecision(int precision);
}
