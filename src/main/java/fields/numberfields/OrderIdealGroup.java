package fields.numberfields;

import java.math.BigInteger;
import java.util.Iterator;

import fields.exceptions.InfinityException;
import fields.integers.Integers;
import fields.interfaces.Group;
import fields.interfaces.Ideal;
import fields.numberfields.NumberField.NFE;
import fields.numberfields.NumberFieldOrder.NumberFieldOrderIdeal;

public class OrderIdealGroup implements Group<FractionalOrderIdeal> {
	private NumberFieldOrder order;

	OrderIdealGroup(NumberFieldOrder order) {
		this.order = order;
	}

	@Override
	public Exactness exactness() {
		return order.exactness();
	}

	@Override
	public FractionalOrderIdeal getRandomElement() {
		NumberFieldOrderIdeal numerator = (NumberFieldOrderIdeal) order.getIdeal(order.getRandomElement(),
				order.getRandomElement());
		NumberFieldOrderIdeal denominator = (NumberFieldOrderIdeal) order.getIdeal(order.getRandomElement(),
				order.getRandomElement());
		return new FractionalOrderIdeal(order, numerator, denominator);
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
	public Iterator<FractionalOrderIdeal> iterator() {
		throw new InfinityException();
	}

	@Override
	public FractionalOrderIdeal neutral() {
		return new FractionalOrderIdeal(order, order.getUnitIdeal(), order.getUnitIdeal());
	}

	public FractionalOrderIdeal getEmbedding(Ideal<NFE> ideal) {
		return new FractionalOrderIdeal(order, (NumberFieldOrderIdeal) ideal);
	}

	public FractionalOrderIdeal getPrincipalIdeal(NFE t) {
		return FractionalOrderIdeal.principalIdeal(order, t);
	}

	@Override
	public FractionalOrderIdeal inverse(FractionalOrderIdeal t) {
		return new FractionalOrderIdeal(order, order.getIdeal(order.getInteger(t.getDenominator())), t.getNumerator());
	}

	@Override
	public FractionalOrderIdeal operate(FractionalOrderIdeal t1, FractionalOrderIdeal t2) {
		return new FractionalOrderIdeal(order,
				(NumberFieldOrderIdeal) order.multiply(t1.getNumerator(), t2.getNumerator()),
				Integers.z().multiply(t1.getDenominator(), t2.getDenominator()));
	}

}
