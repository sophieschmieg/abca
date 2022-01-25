package fields.numberfields;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.TreeMap;

import fields.exceptions.InfinityException;
import fields.helper.AbstractElement;
import fields.helper.AbstractModule;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.interfaces.Group;
import fields.interfaces.Ideal;
import fields.interfaces.Module;
import fields.numberfields.IdealClassGroup.IdealClass;
import fields.numberfields.NumberField.NFE;
import fields.numberfields.NumberFieldIntegers.NumberFieldIdeal;
import fields.vectors.GenericPIDModule;
import fields.vectors.GenericPIDModule.Mod;
import fields.vectors.Vector;

public class IdealClassGroup extends AbstractModule<IntE, IdealClass>
		implements Group<IdealClass>, Module<IntE, IdealClass> {
	private NumberFieldIntegers order;
	private IdealGroup idealGroup;
	private Map<NumberFieldIdeal, IdealClass> reducedMapping;
	private List<IdealClass> elements;
	private List<IdealClass> generators;
	private List<List<IntE>> inGenerators;
	private Map<Vector<IntE>, Integer> generatedIndex;
	private List<List<IntE>> syzygies;
	private int neutralIndex;
	private GenericPIDModule<IntE, Vector<IntE>> asGenericModule;

	public static class IdealClass extends AbstractElement<IdealClass> {
		private NumberFieldIdeal representative;
		private List<NumberFieldIdeal> alternativeRepresentatives;
		private int index;

		private IdealClass(int index, NumberFieldIdeal representative,
				List<NumberFieldIdeal> alternativeRepresentatives) {
			this.index = index;
			this.representative = representative;
			this.alternativeRepresentatives = alternativeRepresentatives;
		}

		@Override
		public int compareTo(IdealClass o) {
			return index - o.index;
		}

		@Override
		public String toString() {
			return "[" + representative + "]";
		}

		public NumberFieldIdeal representative() {
			return representative;
		}

		List<NumberFieldIdeal> alternativeRepresentatives() {
			return alternativeRepresentatives;
		}

		public boolean isPrincipal() {
			return representative.isPrincipal();
		}
	}

	IdealClassGroup(NumberFieldIntegers order) {
		this.order = order;
		this.idealGroup = order.numberField().idealGroup();
		this.reducedMapping = new TreeMap<>();
		this.elements = new ArrayList<>();
		getEmbedding(idealGroup.neutral());
		this.neutralIndex = 0;
	}

	public IdealClass getEmbedding(FractionalIdeal ideal) {
		return getEmbedding(ideal.clearDenominator().getFirst());
	}

	public IdealClass getEmbedding(Ideal<NFE> ideal) {
		NumberFieldIdeal representative = order.reducedRepresentative(ideal).clearDenominator().getFirst();
		if (!reducedMapping.containsKey(representative)) {
			List<NumberFieldIdeal> representatives = order.idealEquivalenceClass(ideal);
			IdealClass idealClass = new IdealClass(elements.size(), representative, representatives);
			elements.add(idealClass);
			for (NumberFieldIdeal alternate : representatives) {
				reducedMapping.put(alternate, idealClass);
			}
		}
		return reducedMapping.get(representative);

	}

	private void computeAllElements() {
		if (generators != null) {
			return;
		}
		Integers z = Integers.z();
		this.generators = new ArrayList<>();
		this.inGenerators = new ArrayList<>();
		this.syzygies = new ArrayList<>();
//		addElement(order.getUnitIdeal());
		this.inGenerators.add(new ArrayList<>());
		IntE minkowskiLimit = order.numberField().minkowskiBound().roundDown();
		for (IntE prime : z.setOfPrimes()) {
			if (prime.compareTo(minkowskiLimit) > 0) {
				break;
			}
			for (NumberFieldIdeal ideal : order.idealsOver(z.getIdeal(prime))) {
				NumberFieldIdeal generatorPower = ideal;
				int power = 1;
				List<IdealClass> elementsCopy = new ArrayList<>();
				elementsCopy.addAll(this.elements);
				boolean generator = false;
				while (!generatorPower.isPrincipal()) {
					AddElementResult addElement = addElement(generatorPower);
					List<IntE> inGenerator = new ArrayList<>();
					int zeroes = generators.size() - 1;
					if (power == 1) {
						zeroes++;
					}
					for (int i = 0; i < zeroes; i++) {
						inGenerator.add(z.zero());
					}
					inGenerator.add(z.getInteger(power));
					if (addElement.newElementFound) {
						System.out.println(generatorPower);
						if (power == 1) {
							generator = true;
							this.generators.add(addElement.element);
							for (List<IntE> g : inGenerators) {
								g.add(z.zero());
							}
							for (List<IntE> syzygy : syzygies) {
								syzygy.add(z.zero());
							}
						}
						inGenerators.add(inGenerator);
						for (IdealClass element : elementsCopy) {
							if (element.representative.equals(order.getUnitIdeal())) {
								continue;
							}
							List<IntE> factorInGenerator = inGenerators.get(element.index);
							AddElementResult multiplied = addElement(
									order.multiply(element.representative, generatorPower));
							if (multiplied.newElementFound) {
								List<IntE> multipliedInGenerator = new ArrayList<>();
								for (int i = 0; i < generators.size(); i++) {
									multipliedInGenerator.add(z.add(factorInGenerator.get(i), inGenerator.get(i)));
								}
								inGenerators.add(multipliedInGenerator);
							}
						}
					} else {
						if (power == 1) {
							break;
						}
						List<IntE> representativeInGenerator = inGenerators.get(addElement.index);
						List<IntE> syzygy = new ArrayList<>();
						for (int i = 0; i < generators.size(); i++) {
							syzygy.add(z.subtract(representativeInGenerator.get(i), inGenerator.get(i)));
						}
						syzygies.add(syzygy);

					}
					generatorPower = order.multiply(ideal, generatorPower);
					power++;
				}
				if (generator) {
					List<IntE> syzygy = new ArrayList<>();
					for (int i = 0; i < generators.size() - 1; i++) {
						syzygy.add(z.zero());
					}
					syzygy.add(z.getInteger(power));
					syzygies.add(syzygy);
				}
			}
		}
		if (!generators.isEmpty()) {
			this.generatedIndex = new TreeMap<>();
			this.asGenericModule = GenericPIDModule.fromSyzygies(Integers.z(), generators.size(), syzygies);
			for (int i = 0; i < elements.size(); i++) {
				generatedIndex.put(
						asGenericModule.asVector(asGenericModule.fromVector(new Vector<>(inGenerators.get(i)))), i);
			}
		}
	}

	@Override
	public String toString() {
		return "Cl(" + order.numberField() + ")";
	}

	@Override
	public Integers getRing() {
		return Integers.z();
	}

	private class AddElementResult {
		private boolean newElementFound;
		private IdealClass element;
		private int index;
	}

	private AddElementResult addElement(NumberFieldIdeal ideal) {
		NumberFieldIdeal representative = order.reducedRepresentative(ideal).clearDenominator().getFirst();
		AddElementResult result = new AddElementResult();
		if (reducedMapping.containsKey(representative)) {
			result.newElementFound = false;
			result.element = reducedMapping.get(representative);
			result.index = result.element.index;
			return result;
		}
//		for (int i = 0; i < elements.size(); i++) {
//			IdealClass element = elements.get(i);
//			if (new FractionalIdeal(order, element.representative, representative).isPrincipal()) {
//				result.newElementFound = false;
//				result.element = element;
//				result.index = i;
//				return result;
//			}
//		}
		result.newElementFound = true;
		result.element = getEmbedding(ideal);// new IdealClass(result.index, representative);
		result.index = result.element.index;
//		elements.add(result.element);
		return result;
	}

	@Override
	public Exactness exactness() {
		return order.exactness();
	}

	@Override
	public IdealClass getRandomElement() {
		return elements.get(new Random().nextInt(elements.size()));
	}

	@Override
	public boolean isFinite() {
		return true;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		computeAllElements();
		return BigInteger.valueOf(elements.size());
	}

	@Override
	public Iterator<IdealClass> iterator() {
		computeAllElements();
		return elements.iterator();
	}

//	public IdealClass getEmbedding(FractionalIdeal ideal) {
//		for (IdealClass element : elements) {
//			if (idealGroup.operate(ideal, idealGroup.inverse(idealGroup.getEmbedding(element.representative)))
//					.isPrincipal()) {
//				return element;
//			}
//		}
//		throw new ArithmeticException("Ideal Class not found!");
//	}
//
//	public IdealClass getEmbedding(NumberFieldIdeal ideal) {
//		return getEmbedding(idealGroup.getEmbedding(ideal));
//	}

	@Override
	public IdealClass neutral() {
		return elements.get(neutralIndex);
	}

	@Override
	public IdealClass zero() {
		return neutral();
	}

	@Override
	public IdealClass inverse(IdealClass t) {
		return getEmbedding(idealGroup.inverse(idealGroup.getEmbedding(t.representative)));
//		if (elements.size() == 1) {
//			return elements.get(0);
//		}
//		return fromIntModuleElement(asGenericModule.negative(asIntModuleElement(t)));
	}

	@Override
	public IdealClass negative(IdealClass s) {
		return inverse(s);
	}

	@Override
	public IdealClass operate(IdealClass t1, IdealClass t2) {
		if (elements.size() == 1) {
			return elements.get(0);
		}
		return getEmbedding(order.multiply(t1.representative, t2.representative));
		// return fromIntModuleElement(asGenericModule.add(asIntModuleElement(t1),
		// asIntModuleElement(t2)));
	}

	@Override
	public IdealClass add(IdealClass s1, IdealClass s2) {
		return operate(s1, s2);
	}

	@Override
	public IdealClass scalarMultiply(IntE t, IdealClass s) {
		return power(t.getValue(), s);
	}

	@Override
	public boolean isFree() {
		return false;
	}

	@Override
	public Ideal<IntE> annihilator() {
		if (elements.size() == 1) {
			return Integers.z().getUnitIdeal();
		}
		return asGenericModule.annihilator();
	}

	@Override
	public boolean isLinearIndependent(List<IdealClass> s) {
		return s.isEmpty();
	}

	private Mod<Vector<IntE>> asIntModuleElement(IdealClass s) {
		return asGenericModule.reduce(asVector(s));
	}

	private IdealClass fromIntModuleElement(Mod<Vector<IntE>> s) {
		return fromVector(asGenericModule.lift(s));
	}

	private List<Mod<Vector<IntE>>> asIntModuleElementList(List<IdealClass> s) {
		List<Mod<Vector<IntE>>> result = new ArrayList<>();
		for (IdealClass c : s) {
			result.add(asIntModuleElement(c));
		}
		return result;
	}

	@Override
	public boolean isGeneratingModule(List<IdealClass> s) {
		if (elements.size() == 1) {
			return true;
		}
		return asGenericModule.isGeneratingModule(asIntModuleElementList(s));
	}

	@Override
	public List<List<IntE>> nonTrivialCombinations(List<IdealClass> s) {
		if (elements.size() == 1) {
			List<List<IntE>> result = new ArrayList<>();
			for (int i = 0; i < s.size(); i++) {
				List<IntE> row = new ArrayList<>();
				for (int j = 0; j < i; j++) {
					row.add(Integers.z().zero());
				}
				row.add(Integers.z().one());
				for (int j = i + 1; j < s.size(); j++) {
					row.add(Integers.z().zero());
				}
				result.add(row);
			}
			return result;
		}
		return asGenericModule.nonTrivialCombinations(asIntModuleElementList(s));
	}

	@Override
	public List<IdealClass> getModuleGenerators() {
		return generators;
	}

	public List<IdealClass> getDiagonalModuleGenerators() {
		if (elements.size() == 1) {
			return Collections.emptyList();
		}
		List<IdealClass> result = new ArrayList<>();
		for (Mod<Vector<IntE>> diagonal : asGenericModule.diagonalBasis()) {
			result.add(fromIntModuleElement(diagonal));
		}
		return result;
	}

	public List<IntE> getDiagonalRanks() {
		if (elements.size() == 1) {
			return Collections.emptyList();
		}
		return asGenericModule.diagonalRanks();
	}

	@Override
	public List<List<IntE>> getModuleGeneratorRelations() {
		return syzygies;
	}

	@Override
	public Vector<IntE> asVector(IdealClass s) {
		return new Vector<>(inGenerators.get(s.index));
	}

	@Override
	public IdealClass fromVector(Vector<IntE> vector) {
		if (elements.size() == 1) {
			return elements.get(0);
		}
		if (!generatedIndex.containsKey(vector)) {
			vector = asGenericModule.asVector(asGenericModule.fromVector(vector));
		}
		return elements.get(generatedIndex.get(vector));
	}
}
