package fields.interfaces;

import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;

public interface ValueField<T extends Element<T>> extends Field<T> {
	boolean isComplete();
	public interface AbsoluteValue<T> {
		public Real value(T t);
	}

	public Real value(T t);

	default public AbsoluteValue<T> value() {
		return new AbsoluteValue<>() {

			@Override
			public Real value(T t) {
				return ValueField.this.value(t);
			}
		};
	}

	public Reals getReals();
}
