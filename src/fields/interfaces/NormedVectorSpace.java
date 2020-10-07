package fields.interfaces;

public interface NormedVectorSpace<T extends Element<T>, S extends Element<S>> extends VectorSpace<T, S> {
	double norm(S s);
	ValueField<T> getValueField();
}
