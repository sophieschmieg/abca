package fields.helper;

import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import fields.interfaces.AlgebraicExtensionElement;
import fields.interfaces.Element;
import fields.interfaces.FieldExtension;
import fields.interfaces.MathMap;
import fields.interfaces.UnivariatePolynomialRing;

public class FieldAutomorphism<T extends Element<T>, S extends AlgebraicExtensionElement<T, S>, Ext extends FieldExtension<T, S, Ext>>
		extends AbstractElement<FieldAutomorphism<T, S, Ext>> implements MathMap<S, S> {
	private int[] images;
	private FieldExtension<T, S, Ext> extension;

	public FieldAutomorphism(FieldExtension<T, S, Ext> extension, int[] images) {
		if (images.length != extension.degree()) {
			throw new ArithmeticException("Dimensions not aligned");
		}
		this.extension = extension;
		this.images = Arrays.copyOf(images, images.length);
		Set<Integer> testSet = new TreeSet<>();
		for (int i = 0; i < this.images.length; i++) {
			if (i < 0 || i >= this.images.length || testSet.contains(this.images[i])) {
				throw new ArithmeticException("Not a permutation");
			}
			testSet.add(this.images[i]);
		}
	}

	public int[] asArray() {
		return images;
	}

	@Override
	public int compareTo(FieldAutomorphism<T, S, Ext> o) {
		for (int i = 0; i < images.length; i++) {
			int cmp = images[i] - o.images[i];
			if (cmp != 0) {
				return cmp;
			}
		}
		return 0;
	}

	@SuppressWarnings("unchecked")
	@Override
	public S evaluate(S t) {
		List<S> baseImages = extension.conjugates(extension.alpha());
		S alphaImage = baseImages.get(images[0]);
		UnivariatePolynomialRing<S> polynomials = extension.getUnivariatePolynomialRing();
		return polynomials.evaluate(polynomials.getEmbedding(t.asPolynomial(), extension.getEmbeddingMap()), alphaImage);
	}
}
