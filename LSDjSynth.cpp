//// LSDjSynth.cpp
//// HertzDevil 2015
//// MIT License.

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <string> // C++ stoi
#include "LSDjSynth.hpp"

using namespace std;

// code begins here

#define max(x, y) ((x) > (y) ? (x) : (y))
#define min(x, y) ((x) < (y) ? (x) : (y))

void abort(err_msg_index id)
{
	const char* err = ERR_MSG[id];
	fprintf(stderr, err != nullptr ? err : "Unknown error");
	putc('\n', stderr);
	exit(EXIT_FAILURE);
}

bool readArgNum(unsigned char* dest, char* str)
{
	try {
		size_t len = 0;
		int v = stoi(str, &len, 0);
		if ((v & ~0xFF) || len < strlen(str))
			return false;
		*dest = static_cast<unsigned char>(v);
	}
	catch (...) {
		return false;
	}
	return true;
}

bool readWave(waveParam* w, char** p)
{
	for (int i = 0; i < 4; i++) {
		if (!readArgNum(&w->volume + i, *p++))
			return false;
	}
	if (w->phase & 0x20)
		return false;
	return true;
}
/*
bool readInst(instParam* inst, char** p)
{
	for (int i = 1; i < 4; i++) {
		unsigned char* dest = reinterpret_cast<unsigned char*>(&inst->play) + i;
		if (!readArgNum(dest, *p++))
			return false;
		if (*dest & 0xF0)
			return false;
		if (i == 3 && !*dest)
			return false;
	}
	return true;
}
*/
void interpolate(unsigned char* dest, size_t len)
{
	bool rev = false;
	unsigned char a = dest[0], b = dest[len - 1];
	if (a == b) {
		memset(dest, a, len);
		return;
	}
	if (rev = a > b) {
		unsigned char c = a; a = b; b = c;
		dest += len - 1;
	}
	if (b != 0xFF) b++;
	for (size_t i = 0; i < len; i++)
		(rev ? *dest-- : *dest++) = a + (i * (b - a) >> 4);
}

const char* getMainWave(wave_t shape, filter_t filter, unsigned char cutoff)
{
	switch (shape) {
	case SAWTOOTH: case SQUARE: case TRIANGLE: break;
	default: abort(E_UNKNOWN_WAVE);
	}
	switch (filter) {
	case LOWP:
		switch (shape) {
		case SAWTOOTH: return SAWTOOTH_LP[cutoff >> 4];
		case SQUARE:   return SQUARE_LP[cutoff >> 5];
		case TRIANGLE: return TRIANGLE_LP[cutoff >> 5];
		}
	case HIGHP:
		switch (shape) {
		case SAWTOOTH: return SAWTOOTH_HP[cutoff >> 4];
		case SQUARE:   return SQUARE_HP[cutoff >> 5];
		case TRIANGLE: return TRIANGLE_HP[cutoff >> 5];
		}
	case BANDP:
		return nullptr;
	case ALLP:
		switch (shape) {
		case SAWTOOTH: return reinterpret_cast<const char*>(SAWTOOTH_LP) + (sizeof(SAWTOOTH_LP) - WAVE_LENGTH) / 2;
		case SQUARE:   return reinterpret_cast<const char*>(SQUARE_LP) + (sizeof(SQUARE_LP) - WAVE_LENGTH) / 2;
		case TRIANGLE: return reinterpret_cast<const char*>(TRIANGLE_LP) + (sizeof(TRIANGLE_LP) - WAVE_LENGTH) / 2;
		}
	}
	abort(E_UNKNOWN_FILTER);
	return nullptr;
}

const char* getPartial(wave_t shape, unsigned char cutoff)
{
	switch (shape) {
	case SAWTOOTH: return PARTIAL[cutoff >> 4];
	case SQUARE:   return PARTIAL[cutoff >> 5 << 1];
	case TRIANGLE: return PARTIAL_TRI[min(4U, cutoff >> 5)];
	}
	abort(E_UNKNOWN_WAVE);
	return nullptr;
}

const char* getSidePartial(wave_t shape, unsigned char cutoff)
{
	switch (shape) {
	case SAWTOOTH: return cutoff < 0x10 ? nullptr : PARTIAL[(cutoff >> 4) - 1];
	case SQUARE:   return cutoff < 0x20 ? nullptr : PARTIAL[(cutoff >> 5 << 1) - 1];
	case TRIANGLE: return cutoff < 0x20 ? nullptr : PARTIAL_TRI[min(4U, (cutoff >> 5) - 1)];
	}
	abort(E_UNKNOWN_WAVE);
	return nullptr;
}

unsigned char getPartialVelocity(wave_t shape, unsigned char cutoff)
{
	switch (shape) {
	case SAWTOOTH: return cutoff & 0x0F;
	case SQUARE:
	case TRIANGLE: return (cutoff >> 1) & 0x0F;
	}
	abort(E_UNKNOWN_WAVE);
	return 0U;
}

static bool Intp = false;

bool makeWaves(int (table[])[WAVE_LENGTH], synthParam* synth)
{
	const unsigned int LEN = Intp ? WAVE_COUNT : 1;

	unsigned char* intp_volume = new unsigned char[LEN];
	unsigned char* intp_cutoff = new unsigned char[LEN];
	unsigned char* intp_phase  = new unsigned char[LEN];
	unsigned char* intp_vshift = new unsigned char[LEN];
	
	intp_volume[0] = synth->start.volume;
	intp_cutoff[0] = synth->start.cutoff;
	intp_phase[0]  = synth->start.phase;
	intp_vshift[0] = synth->start.vshift;
	intp_volume[LEN - 1] = synth->end.volume;
	intp_cutoff[LEN - 1] = synth->end.cutoff;
	intp_phase[LEN - 1]  = synth->end.phase;
	intp_vshift[LEN - 1] = synth->end.vshift;
	
	interpolate(intp_volume, LEN);
	interpolate(intp_cutoff, LEN);
	interpolate(intp_phase, LEN);
	interpolate(intp_vshift, LEN);

	int qCoeff, qsCoeff;
	switch (synth->filter) {
	case LOWP:  qCoeff = synth->q + 1;            qsCoeff = synth->q; break;
	case HIGHP: qCoeff = synth->q    ;            qsCoeff = synth->q + 1; break;
	case BANDP: qCoeff = synth->q ? synth->q : 1; qsCoeff = qCoeff; break;
	case ALLP:  qCoeff = synth->q;                qsCoeff = qCoeff; break;
	default: abort(E_UNKNOWN_FILTER);
	}

	int temp[WAVE_COUNT][WAVE_LENGTH] = {};
	for (size_t i = 0; i < LEN; i++) {
		const char* mainWave = getMainWave(synth->wave, synth->filter, intp_cutoff[i]);
		const char* partial = getPartial(synth->wave, intp_cutoff[i]);
		const char* side = getSidePartial(synth->wave, intp_cutoff[i]);
		const unsigned char vel = getPartialVelocity(synth->wave, intp_cutoff[i]);
		for (size_t j = 0; j < WAVE_LENGTH; j++) {
			int val = 0;
			if (mainWave != nullptr) val += mainWave[j];
			if (partial  != nullptr) val += partial[j]  * vel          * qCoeff  >> 4;
			if (side     != nullptr) val += side[j]     * (0x0F - vel) * qsCoeff >> 4;
			val = val * intp_volume[i] + 0x800;
			if (val >= 0x4000) val = 0;
			temp[i][j] = ((synth->dist == CLIP ? max(0, min(0xFFF, val)) : val) >> 4) + intp_vshift[i];
		}
	}
	
	for (size_t i = 0; i < LEN; i++) {
		unsigned int wl = 0x20 - intp_phase[i];
		for (size_t j = 0; j < WAVE_LENGTH; j++) {
			size_t index = 0;
			switch (synth->phase) {
			case NORMAL: index = 0x20 * min(j, wl - 1) / wl; break;
			case RESYNC: index = 0x20 * (j % wl) / wl; break;
			case RESYN2: index = j % wl; break;
			default: abort(E_UNKNOWN_PHASE);
			}
			table[i][j] = temp[i][index];
		}
	}

	return true;
}

bool reduceWaves(int (table[])[WAVE_LENGTH], synthParam* synth, waveTable out)
{
	for (size_t i = 0; i < (Intp ? WAVE_COUNT : 1U); i++) {
		for (size_t j = 0; j < WAVE_LENGTH; j++) {
			out[i][j] = table[i][j] >> 4 & 0x0F;
			if (j & 0x01)
				out[i][j - 1] = (out[i][j - 1] + (table[i][j] >> 8 & 0x0F)) & 0x0F;
		}
	}

	return true;
}

sequence_t* makeSequence(instParam* inst)
{
	sequence_t* seq = new sequence_t();

	switch (inst->play) {
	case ONCE:
		seq->length = inst->speed * inst->length + 1;
		seq->loop = -1;
		seq->data[inst->speed * inst->length] = inst->length;
		break;
	case LOOP:
		seq->length = inst->speed * (inst->length + 1);
		seq->loop = inst->speed * (inst->length - inst->repeat);
		break;
	case PINGPONG:
		if ((seq->length = inst->speed * (inst->length + max(1, inst->repeat))) > 252)
			abort(E_FTI_SEQ_TOO_LONG);
		seq->loop = inst->speed * (inst->length - inst->repeat);
		for (unsigned char i = 1; i < inst->repeat; i++)
			memset(seq->data + (i + inst->length) * inst->speed, i ? (0x0F * (inst->length - i) / inst->length) : 0, inst->speed);
		break;
	case MANUAL:
		seq->length = 0;
		seq->loop = -1;
		break;
	default:
		abort(E_PLAY);
	}
	if (inst->play != MANUAL)
		for (unsigned char i = 0; i <= inst->length; i++)
			memset(seq->data + i * inst->speed, i ? (0x0F * i / inst->length) : 0, inst->speed);

	return seq;
}

FILE* createFileExt(const char* fname, const char* ext, const char* mode)
{
	size_t len = strlen(fname);
	const size_t tot = len + strlen(ext) + 2;
	char* buf = new char[tot]();
	if (strcpy_s(buf, tot, fname)) abort(E_FILE);
	if (strcat_s(buf, tot, ".")) abort(E_FILE);
	if (strcat_s(buf, tot, ext)) abort(E_FILE);
	FILE* f = nullptr;
	if (fopen_s(&f, buf, mode)) abort(E_FILE);
	delete[] buf;
	if (!f) abort(E_FILE);
	return f;
}

#define CHECK { if (ferror(f)) { fclose(f); return false; } }

bool createBIN(waveTable table, bool intp, char* fname)
{
	FILE* f = createFileExt(fname, "bin", "wb");

	for (size_t i = 0; i < (intp ? WAVE_COUNT : 1U); i++)
		for (size_t j = 0; j < WAVE_LENGTH; j += 2) {
			CHECK fputc((table[i][j + 1] & 0x0F) | (table[i][j] << 4 & 0xF0), f);
		}
	
	fclose(f);
	return true;
}

bool createFTI(waveTable table, bool intp, char* fname, instParam* inst)
{
	FILE* f = createFileExt(fname, "fti", "wb");

	sequence_t* seq = makeSequence(inst);
	int x = 0;

	static const char FTI_HEADER[] = {
		'F', 'T', 'I', '2', '.', '4',		// identifier
		0x05,								// chip type (N163)
		0x09, 0x00, 0x00, 0x00,				// instrument name length
		'L', 'S', 'D', 'j', 'S', 'y', 'n', 't', 'h',
		0x05,								// number of sequences
	};
	CHECK fwrite(FTI_HEADER, sizeof(char), sizeof(FTI_HEADER), f);

	if (inst->play == ONCE) {
		CHECK fputc(1, f);
		int length = inst->speed * (inst->length + 1) + 1;
		CHECK fwrite(reinterpret_cast<char*>(&length), sizeof(int), 1, f);			// length
		CHECK fwrite(reinterpret_cast<char*>(&(x = -1)), sizeof(int), 1, f);		// loop
		CHECK fwrite(reinterpret_cast<char*>(&(x = -1)), sizeof(int), 1, f);		// release
		CHECK fwrite(reinterpret_cast<char*>(&(x = 0)), sizeof(int), 1, f);			// setting
		for (int i = 0; i < length - 1; i++) {
			CHECK fputc(0x0F, f);			// volume
		}
	}
	CHECK fputc(0, f);						// volume seq last term, or volume seq disable

	CHECK fputc(0, f);						// arp
	CHECK fputc(0, f);						// pitch
	CHECK fputc(0, f);						// hi-pitch
	
	if (intp && inst->play != MANUAL) {
		CHECK fputc(1, f);					// duty
		CHECK fwrite(reinterpret_cast<char*>(&seq->length), sizeof(int), 1, f);		// length
		CHECK fwrite(reinterpret_cast<char*>(&(seq->loop)), sizeof(int), 1, f);		// loop
		CHECK fwrite(reinterpret_cast<char*>(&(x = -1)), sizeof(int), 1, f);		// release
		CHECK fwrite(reinterpret_cast<char*>(&(x = 0)), sizeof(int), 1, f);			// setting
		CHECK fwrite(seq->data, sizeof(char), seq->length, f);
	}
	else {
		CHECK fputc(0, f);					// duty
	}
	
	CHECK fwrite(reinterpret_cast<char*>(&(x = WAVE_LENGTH)), sizeof(int), 1, f);
	CHECK fwrite(reinterpret_cast<char*>(&(x = 0)), sizeof(int), 1, f);
	CHECK fwrite(reinterpret_cast<char*>(&(x = (intp ? WAVE_COUNT : 1))), sizeof(int), 1, f);
	for (size_t i = 0; i < (intp ? WAVE_COUNT : 1U); i++)
		for (size_t j = 0; j < WAVE_LENGTH; j++) {
			CHECK fputc(table[i][j], f);
		}

	fclose(f);
	return true;
}

bool createXPMCK(waveTable table, bool intp, char* fname, instParam* inst)
{
	FILE* f = createFileExt(fname, "mml", "w");
	
	for (unsigned char i = 0; i <= (intp ? inst->length : 0); i++) {
		CHECK fprintf(f, "@WT%d = { ", i);
		const size_t index = i ? (0x0F * i / inst->length) : 0;
		for (size_t j = 0; j < WAVE_LENGTH; j++) {
			if (ferror(f)) return false;
			CHECK fprintf(f, "%d ", table[index][j]);
		}
		CHECK fprintf(f, "}\n");
	}

	if (!intp || inst->play == MANUAL) {
		fclose(f);
		return true;
	}
	CHECK fprintf(f, "@WTM0 = { ");
	for (int i = 0; i <= inst->length; i++) {
		if (i == inst->length - inst->repeat && inst->play != ONCE) {
			CHECK fprintf(f, "| ");
		}
		CHECK fprintf(f, "WT%d %d ", i, inst->speed);
	}
	if (inst->play == PINGPONG)
		for (int i = 1; i < inst->repeat; i++) {
			CHECK fprintf(f, "WT%d %d ", inst->length - i, inst->speed);
		}
	CHECK fprintf(f, "}\n");
	
	fclose(f);
	return true;
}

#undef CHECK

void printWaves(waveTable table, bool intp)
{
	for (size_t i = 0; i < (intp ? WAVE_COUNT : 1U); i++) {
		printf("WAVE %X: ", i);
		for (size_t j = 0; j < WAVE_LENGTH; j++)
			putchar(table[i][j] + (table[i][j] > 9 ? '7' : '0'));
		putchar('\n');
	}
}

int main(int argc, char** argv)
{
	for (int i = 1; i < argc; i++) {
		if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-?")) {
			printf(HELP, argv[0]);
			exit(EXIT_SUCCESS);
		}
		else if (!strcmp(argv[i], "--version") || !strcmp(argv[i], "-v")) {
			printf(VERSION);
			exit(EXIT_SUCCESS);
		}
	}

	if (argc < 5) {
		fprintf(stderr, HELP, argv[0]);
		exit(EXIT_SUCCESS);
	}

	synthParam synth = DEFAULT_SYNTH;
	instParam waveInst = DEFAULT_INST;
	int exportMode = EXPORT_NONE;
	char* fname = "lsdjsynth_out";
	if (!readWave(&synth.start, argv + 1))
		abort(E_SYNTH_PARAM);
	synth.end = synth.start;
	
	const char* ARG[] = {
		"sawtooth", "square", "triangle",
		"lowp", "highp", "bandp", "allp",
		"q",
		"clip", "wrap",
		"normal", "resync", "resyn2",
	};

	for (int i = 5; i < argc; i++) {
		if (!strcmp(argv[i], "--sawtooth") || !strcmp(argv[i], "-saw"))
			synth.wave = SAWTOOTH;
		else if (!strcmp(argv[i], "--square") || !strcmp(argv[i], "-sq"))
			synth.wave = SQUARE;
		else if (!strcmp(argv[i], "--triangle") || !strcmp(argv[i], "-tri"))
			synth.wave = TRIANGLE;

		else if (!strcmp(argv[i], "--lowp") || !strcmp(argv[i], "-lp"))
			synth.filter = LOWP;
		else if (!strcmp(argv[i], "--highp") || !strcmp(argv[i], "-hp"))
			synth.filter = HIGHP;
		else if (!strcmp(argv[i], "--bandp") || !strcmp(argv[i], "-bp"))
			synth.filter = BANDP;
		else if (!strcmp(argv[i], "--allp") || !strcmp(argv[i], "-all"))
			synth.filter = ALLP;

		else if (!strcmp(argv[i], "--q") || !strcmp(argv[i], "-q")) {
			if (++i >= argc)
				abort(E_Q_PARAM);
			if (!readArgNum(&synth.q, argv[i]))
				abort(E_Q_PARAM);
			if (synth.q & 0xF0)
				abort(E_Q_PARAM);
		}

		else if (!strcmp(argv[i], "--clip"))
			synth.dist = CLIP;
		else if (!strcmp(argv[i], "--wrap") || !strcmp(argv[i], "-w"))
			synth.dist = WRAP;

		else if (!strcmp(argv[i], "--normal"))
			synth.phase = NORMAL;
		else if (!strcmp(argv[i], "--resync") || !strcmp(argv[i], "-r"))
			synth.phase = RESYNC;
		else if (!strcmp(argv[i], "--resyn2") || !strcmp(argv[i], "-r2"))
			synth.phase = RESYN2;

		else if (!strcmp(argv[i], "--once") || !strcmp(argv[i], "-o"))
			waveInst.play = ONCE;
		else if (!strcmp(argv[i], "--loop") || !strcmp(argv[i], "-l"))
			waveInst.play = LOOP;
		else if (!strcmp(argv[i], "--pingpong") || !strcmp(argv[i], "-p"))
			waveInst.play = PINGPONG;
		else if (!strcmp(argv[i], "--manual") || !strcmp(argv[i], "-m"))
			waveInst.play = MANUAL;

		else if (!strcmp(argv[i], "--length") || !strcmp(argv[i], "-len")) {
			if (++i >= argc)
				abort(E_INST_PARAM);
			if (!readArgNum(&waveInst.length, argv[i]))
				abort(E_INST_PARAM);
			if (waveInst.length & 0xF0)
				abort(E_INST_PARAM);
		}

		else if (!strcmp(argv[i], "--repeat") || !strcmp(argv[i], "-rep")) {
			if (++i >= argc)
				abort(E_INST_PARAM);
			if (!readArgNum(&waveInst.repeat, argv[i]))
				abort(E_INST_PARAM);
			if (waveInst.repeat & 0xF0)
				abort(E_INST_PARAM);
		}

		else if (!strcmp(argv[i], "--speed") || !strcmp(argv[i], "-spd")) {
			if (++i >= argc)
				abort(E_INST_PARAM);
			if (!readArgNum(&waveInst.speed, argv[i]))
				abort(E_INST_PARAM);
			if ((waveInst.speed & 0xF0) || !waveInst.speed)
				abort(E_INST_PARAM);
		}

		else if (!strcmp(argv[i], "--end") || !strcmp(argv[i], "-e")) {
			if (++i + 4 > argc)
				abort(E_SYNTH_PARAM);
			if (!readWave(&synth.end, (argv + i)))
				abort(E_SYNTH_PARAM);
			i += 3;
			Intp = true;
		}
		
		else if (!strcmp(argv[i], "--bin") || !strcmp(argv[i], "-b"))
			exportMode |= EXPORT_BIN;
		else if (!strcmp(argv[i], "--fti") || !strcmp(argv[i], "-f"))
			exportMode |= EXPORT_FTI;
		else if (!strcmp(argv[i], "--xpmck") || !strcmp(argv[i], "-x"))
			exportMode |= EXPORT_XPMCK;
		
		else if (!strcmp(argv[i], "--fname") || !strcmp(argv[i], "-n"))
			fname = argv[++i];

		else
			abort(E_UNKNOWN_OPTION);
	}

	waveInst.repeat = min(waveInst.repeat, waveInst.length);

	int temp[WAVE_COUNT][WAVE_LENGTH] = {};
	makeWaves(temp, &synth);
	unsigned char raw[WAVE_COUNT][WAVE_LENGTH] = {};
	reduceWaves(temp, &synth, raw);
	if (exportMode & EXPORT_BIN)
		if (!createBIN(raw, Intp, fname))
			abort(E_EXPORT);
	if (exportMode & EXPORT_FTI)
		if (!createFTI(raw, Intp, fname, &waveInst))
			abort(E_EXPORT);
	if (exportMode & EXPORT_XPMCK)
		if (!createXPMCK(raw, Intp, fname, &waveInst))
			abort(E_EXPORT);
	//if (!exportMode)
		printWaves(raw, Intp);

    exit(EXIT_SUCCESS); // return 0;
}
