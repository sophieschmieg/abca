package fields.interfaces;

public interface DivisionAlgebra<T extends Element<T>, S extends Element<S>> extends Algebra<T, S>, DivisionRing<S> {
	Field<T> getField();

	default Ring<T> getRing() {
		return getField();
	}
}
