package cryptography;

import java.util.Arrays;

import fields.helper.AbstractElement;

public class ByteArray extends AbstractElement<ByteArray> {
	private byte[] array;

	public static ByteArray concatenate(byte[] first, byte[] second) {
		byte[] array = new byte[first.length + second.length];
		System.arraycopy(first, 0, array, 0, first.length);
		System.arraycopy(second, 0, array, first.length, second.length);
		return new ByteArray(array);
	}

	public static ByteArray concatenate(ByteArray first, ByteArray second) {
		return concatenate(first.array, second.array);
	}

	public ByteArray(byte[] array) {
		this.array = array;
	}

	public byte[] array() {
		return array;
	}

	public byte byteAt(int index) {
		return array[index];
	}

	public boolean bitAt(int index) {
		return (array[index / 8] & (1 << (index % 8))) != 0;
	}

	@Override
	public int compareTo(ByteArray o) {
		return Arrays.compare(array, o.array);
	}

	@Override
	public String toString() {
		String[] hex = new String[array.length];
		for (int i = 0; i < array.length; i++) {
			int b = array[i] & 0xff;
			int digit1 = b % 16;
			int digit16 = (b / 16) % 16;
			hex[i] = Integer.toHexString(digit16) + Integer.toHexString(digit1);
		}
		return String.join(" ", hex);
	}

	public int length() {
		return array.length;
	}
}
