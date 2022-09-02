package util;

import java.util.function.Function;

import fields.interfaces.Element;
import fields.interfaces.MathMap;

public class FunctionMathMap<T extends Element<T>, S extends Element<S>> implements MathMap<T, S> {
	private Function<T, S> function;

	public FunctionMathMap(Function<T, S> function) {
		this.function = function;
	}

	@Override
	public S evaluate(T t) {
		return function.apply(t);
	}

}
