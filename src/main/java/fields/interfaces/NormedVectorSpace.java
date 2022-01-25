package fields.interfaces;

import fields.floatingpoint.Reals.Real;

public interface NormedVectorSpace<T extends Element<T>, S extends Element<S>> extends VectorSpace<T, S> {
	Real valueNorm(S s);
	Real inverseValueNorm(S s);
	S normedVector(S s);
	ValueField<T> getValueField();
}
