package fields.interfaces;

import fields.floatingpoint.Reals.Real;
import fields.interfaces.ValueField.AbsoluteValue;

public interface ValueRing<T extends Element<T>> extends Ring<T> {
	public Real value(T t);
	public AbsoluteValue<T> value();
}
