package util;

import java.io.IOException;
import java.io.Reader;
import java.util.Optional;

public class PeekableReader extends Reader {
	private final Reader in;
	private boolean closed;
	private StringBuilder buffer;

	public PeekableReader(Reader in) {
		this.in = in;
		this.buffer = new StringBuilder();
		this.closed = false;
	}

	public int peek() throws IOException {
		if (buffer.length() > 0) {
			return buffer.charAt(0);
		}
		int result = in.read();
		if (result < 0) {
			return result;
		}
		buffer.append((char) result);
		return result;
	}

	@Override
	public int read() throws IOException {
		int result;
		if (buffer.length() > 0) {
			result = buffer.charAt(0);
			buffer.delete(0, 1);
			return result;
		}
		return in.read();
	}

	@Override
	public int read(char[] cbuf, int off, int len) throws IOException {
		if (closed) {
			throw new IOException("Closed Reader!");
		}
		int read = 0;
		for (int i = 0; i < Math.min(len, buffer.length()); i++) {
			cbuf[off + i] = buffer.charAt(i);
			read++;
		}
		buffer.delete(0, read);
		if (read < len) {
			int num = in.read(cbuf, off + read, len - read);
			if (num == -1) {
				if (read == -1) {
					return -1;
				}
				return read;
			}
			read += num;
		}
		return read;
	}

	public int peek(char[] cbuf) throws IOException {
		return peek(cbuf, 0, cbuf.length);
	}

	public int peek(char[] cbuf, int off, int len) throws IOException {
		if (buffer.length() < len) {
			char[] data = new char[len - buffer.length()];
			int result = in.read(data, 0, len - buffer.length());
			if (result < 0) {
				if (buffer.length() == 0) {
					return -1;
				}
			} else {
				buffer.append(data);
			}
		}
		buffer.getChars(0, Math.min(len, buffer.length()), cbuf, off);
		return Math.min(len, buffer.length());
	}

	public int peek(int n) throws IOException {
		if (buffer.length() <= n) {
			char[] toRead = new char[n - buffer.length() + 1];
			int read = in.read(toRead, 0, n - buffer.length() + 1);
			if (read < 0) {
				return -1;
			}
			buffer.append(toRead, 0, read);
		}
		return buffer.charAt(n);
	}

	@Override
	public boolean ready() throws IOException {
		return buffer.length() > 0 || in.ready();
	}

	@Override
	public long skip(long n) throws IOException {
		long skipped = Math.min(n, buffer.length());
		buffer.delete(0, (int) skipped);
		if (skipped < n) {
			skipped += in.skip(n - skipped);
		}
		return skipped;
	}

	@Override
	public void close() throws IOException {
		in.close();
		this.closed = true;
	}
}
