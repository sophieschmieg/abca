package fields.helper;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import fields.exceptions.InfinityException;
import fields.interfaces.AlgebraicExtensionElement;
import fields.interfaces.Element;
import fields.interfaces.FieldExtension;
import fields.interfaces.Group;

public class GaloisGroup<T extends Element<T>, S extends AlgebraicExtensionElement<T, S>, Ext extends FieldExtension<T, S, Ext>>
		implements Group<FieldAutomorphism<T, S, Ext>> {
	private FieldExtension<T, S, Ext> fieldExtension;
	private List<FieldAutomorphism<T, S, Ext>> elements;
	private int order;


	public GaloisGroup(FieldExtension<T, S, Ext> fieldExtension, FieldAutomorphism<T, S, Ext> generator) {
		this.fieldExtension = fieldExtension;
		this.order = fieldExtension.degree();
			this.elements = new ArrayList<>();
		FieldAutomorphism<T, S, Ext> element = neutral();
		do  {
			this.elements.add(element);
			element = operate(element, generator);
		} while (!element.equals(neutral()));
	}

	public GaloisGroup(FieldExtension<T, S, Ext> fieldExtension, List<FieldAutomorphism<T, S, Ext>> elements) {
		this.fieldExtension = fieldExtension;
		this.elements = elements;
		this.order = fieldExtension.degree();
	}

	@Override
	public Exactness exactness() {
		return fieldExtension.exactness();
	}

	@Override
	public FieldAutomorphism<T, S, Ext> getRandomElement() {
		Random rng = new Random();
		return elements.get(rng.nextInt(elements.size()));
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
	public Iterator<FieldAutomorphism<T, S, Ext>> iterator() {
		return elements.iterator();
	}

	@Override
	public FieldAutomorphism<T, S, Ext> neutral() {
		int[] images = new int[order];
		for (int i = 0; i < order; i++) {
			images[i] = i;
		}
		return new FieldAutomorphism<>(fieldExtension, images);
	}

	@Override
	public FieldAutomorphism<T, S, Ext> inverse(FieldAutomorphism<T, S, Ext> t) {
		int[] images = new int[order];
		for (int i = 0; i < order; i++) {
			images[t.asArray()[i]] = i;
		}
		return new FieldAutomorphism<>(fieldExtension, images);
	}

	@Override
	public FieldAutomorphism<T, S, Ext> operate(FieldAutomorphism<T, S, Ext> t1, FieldAutomorphism<T, S, Ext> t2) {
		int[] images = new int[order];
		for (int i = 0; i < order; i++) {
			images[i] = t1.asArray()[t2.asArray()[i]];
		}
		return new FieldAutomorphism<>(fieldExtension, images);
	}

}