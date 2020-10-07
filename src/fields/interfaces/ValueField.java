package fields.interfaces;

public interface ValueField<T extends Element<T>> extends Field<T> {
	public double value(T t);
}
