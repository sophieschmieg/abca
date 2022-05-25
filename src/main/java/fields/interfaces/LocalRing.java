package fields.interfaces;

public interface LocalRing<T extends Element<T>, S extends Element<S>, R extends Element<R>>
		extends Ring<T>, DedekindRing<T, S, R> {
	public Ideal<T> maximalIdeal();

	public Field<R> reduction();
	public R reduce(T t);
	public T lift(R t);

	public ModuloMaximalIdealResult<T, R, ?, ?, ?> moduloMaximalIdeal(Ideal<T> ideal);
	public FieldOfFractionsResult<T, S> fieldOfFractions();
}
