package fields.numberfields;

import java.math.BigInteger;
import java.util.Iterator;

import fields.exceptions.InfinityException;
import fields.interfaces.Group;
import fields.numberfields.NumberFieldIntegers.NumberFieldIdeal;

public class IdealGroup implements Group<FractionalIdeal> {
	private NumberFieldIntegers order;

	IdealGroup(NumberFieldIntegers order) {
		this.order = order;
	}

	@Override
	public Exactness exactness() {
		return order.exactness();
	}

	@Override
	public FractionalIdeal getRandomElement() {
		NumberFieldIdeal numerator = (NumberFieldIdeal) order.getIdeal(order.getRandomElement(),
				order.getRandomElement());
		NumberFieldIdeal denominator = (NumberFieldIdeal) order.getIdeal(order.getRandomElement(),
				order.getRandomElement());
		return new FractionalIdeal(order, numerator, denominator);
	}

	@Override
	public boolean isFinite() {
		return false;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		throw new InfinityException();
	}

	@Override
	public Iterator<FractionalIdeal> iterator() {
		throw new InfinityException();
	}

	@Override
	public FractionalIdeal neutral() {
		return new FractionalIdeal(order, order.getUnitIdeal(), order.getUnitIdeal());
	}
	
	public FractionalIdeal getEmbedding(NumberFieldIdeal ideal ) {
		return new FractionalIdeal(order, ideal);
	}

	@Override
	public FractionalIdeal inverse(FractionalIdeal t) {
		return new FractionalIdeal(order, t.getDenominator(), t.getNumerator());
	}

	@Override
	public FractionalIdeal operate(FractionalIdeal t1, FractionalIdeal t2) {
		return new FractionalIdeal(order, order.multiply(t1.getNumerator(), t2.getNumerator()),
				order.multiply(t1.getDenominator(), t2.getDenominator()));
	}

}
