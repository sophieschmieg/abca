package fields.numberfields;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.TreeMap;

import fields.exceptions.InfinityException;
import fields.helper.AbstractElement;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.interfaces.Group;
import fields.numberfields.IdealClassGroup.IdealClassGroupElement;
import fields.numberfields.NumberFieldIntegers.NumberFieldIdeal;

public class IdealClassGroup implements Group<IdealClassGroupElement> {
	private NumberFieldIntegers order;
	private IdealGroup idealGroup;
	private List<IdealClassGroupElement> elements;
	private int neutralIndex;
	private Map<Integer, Map<Integer, Integer>> multiplicationMap;
	private Map<Integer, Integer> inverseMap;

	public static class IdealClassGroupElement extends AbstractElement<IdealClassGroupElement> {
		private NumberFieldIdeal representative;
		private int index;

		private IdealClassGroupElement(int index, NumberFieldIdeal representative) {
			this.index = index;
			this.representative = representative;
		}

		@Override
		public int compareTo(IdealClassGroupElement o) {
			return index - o.index;
		}

		@Override
		public String toString() {
			return "[" + representative + "]";
		}
	}

	IdealClassGroup(NumberFieldIntegers order) {
		Integers z = Integers.z();
		this.order = order;
		this.idealGroup = order.numberField().idealGroup();
		this.elements = new ArrayList<>();
		addElement(order.getUnitIdeal());
		this.neutralIndex = 0;
		this.multiplicationMap = new TreeMap<>();
		this.inverseMap = new TreeMap<>();
		IntE minkowskiLimit = order.numberField().minkowskiBound().roundDown();
		for (IntE prime : z.setOfPrimes()) {
			if (prime.compareTo(minkowskiLimit) > 0) {
				break;
			}
			for (NumberFieldIdeal ideal : order.idealsOver(z.getIdeal(prime))) {
				NumberFieldIdeal generatorPower = ideal;
				int power = 1;
				List<IdealClassGroupElement> elementsCopy = new ArrayList<>();
				elementsCopy.addAll(this.elements);
				while (!generatorPower.isPrincipal()) {
					if (addElement(generatorPower)) {
						for (IdealClassGroupElement element : elementsCopy) {
							if (element.representative.equals(order.getUnitIdeal())) {
								continue;
							}
							addElement(order.multiply(element.representative, generatorPower));
						}
					} else if (power == 1) {
						break;
					}
					generatorPower = order.multiply(ideal, generatorPower);
					power++;
				}
			}
		}
	}

	private boolean addElement(NumberFieldIdeal representative) {
		for (IdealClassGroupElement element : elements) {
			if (new FractionalIdeal(order, element.representative, representative).isPrincipal()) {
				return false;
			}
		}
		int index = elements.size();
		elements.add(new IdealClassGroupElement(index, representative));
		return true;
	}

	@Override
	public Exactness exactness() {
		return order.exactness();
	}

	@Override
	public IdealClassGroupElement getRandomElement() {
		return elements.get(new Random().nextInt(elements.size()));
	}

	@Override
	public boolean isFinite() {
		return true;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		return BigInteger.valueOf(elements.size());
	}

	@Override
	public Iterator<IdealClassGroupElement> iterator() {
		return elements.iterator();
	}

	public IdealClassGroupElement getEmbedding(FractionalIdeal ideal) {
		for (IdealClassGroupElement element : elements) {
			if (idealGroup.operate(ideal, idealGroup.inverse(idealGroup.getEmbedding(element.representative)))
					.isPrincipal()) {
				return element;
			}
		}
		throw new ArithmeticException("Ideal Class not found!");
	}

	public IdealClassGroupElement getEmbedding(NumberFieldIdeal ideal) {
		return getEmbedding(idealGroup.getEmbedding(ideal));
	}

	@Override
	public IdealClassGroupElement neutral() {
		return elements.get(neutralIndex);
	}

	@Override
	public IdealClassGroupElement inverse(IdealClassGroupElement t) {
		if (!inverseMap.containsKey(t.index)) {
			int index = getEmbedding(idealGroup.inverse(idealGroup.getEmbedding(t.representative))).index;
			inverseMap.put(t.index, index);
		}
		return elements.get(inverseMap.get(t.index));
	}

	@Override
	public IdealClassGroupElement operate(IdealClassGroupElement t1, IdealClassGroupElement t2) {
		if (!multiplicationMap.containsKey(t1.index)) {
			multiplicationMap.put(t1.index, new TreeMap<>());
		}
		Map<Integer, Integer> map = multiplicationMap.get(t1.index);
		if (!map.containsKey(t2.index)) {
			int index = getEmbedding(order.multiply(t1.representative, t2.representative)).index;
			map.put(t2.index, index);
		}
		return elements.get(map.get(t2.index));
	}
}
