package fields.interfaces;

public interface MathMap<T extends Element<T>, S extends Element<S>> {
	S evaluate(T t);
}
