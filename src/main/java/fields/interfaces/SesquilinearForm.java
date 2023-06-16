package fields.interfaces;

public interface SesquilinearForm<T extends Element<T>, S extends Element<S>> {
	T evaluate(S t1, S t2);
}
