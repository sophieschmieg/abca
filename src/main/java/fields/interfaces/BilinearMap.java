package fields.interfaces;

public interface BilinearMap<T extends Element<T>, S extends Element<S>> {
	S evaluate(T t1, T t2);
}
