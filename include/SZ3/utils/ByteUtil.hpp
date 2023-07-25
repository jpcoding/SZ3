//
// Created by Kai Zhao on 1/28/21.
//

#ifndef SZ3_BYTEUTIL_HPP
#define SZ3_BYTEUTIL_HPP

#include "SZ3/def.hpp"
#include <cstring>
#include <stdio.h>

namespace SZ {

    typedef union lint16 {
        unsigned short usvalue;
        short svalue;
        unsigned char byte[2];
    } lint16;

    typedef union lint32 {
        int ivalue;
        unsigned int uivalue;
        unsigned char byte[4];
    } lint32;

    typedef union lint64 {
        int64_t lvalue;
        uint64_t ulvalue;
        unsigned char byte[8];
    } lint64;

    typedef union ldouble {
        double value;
        uint64_t lvalue;
        unsigned char byte[8];
    } ldouble;

    typedef union lfloat {
        float value;
        unsigned int ivalue;
        unsigned char byte[4];
        uint16_t int16[2];
    } lfloat;

    inline void symTransform_4bytes(uchar data[4]) {
        unsigned char tmp = data[0];
        data[0] = data[3];
        data[3] = tmp;

        tmp = data[1];
        data[1] = data[2];
        data[2] = tmp;
    }

    inline int16_t bytesToInt16_bigEndian(unsigned char *bytes) {
        int16_t temp = 0;
        int16_t res = 0;

        temp = bytes[0] & 0xff;
        res |= temp;

        res <<= 8;
        temp = bytes[1] & 0xff;
        res |= temp;

        return res;
    }

    inline int32_t bytesToInt32_bigEndian(const unsigned char *bytes) {
        int32_t temp = 0;
        int32_t res = 0;

        res <<= 8;
        temp = bytes[0] & 0xff;
        res |= temp;

        res <<= 8;
        temp = bytes[1] & 0xff;
        res |= temp;

        res <<= 8;
        temp = bytes[2] & 0xff;
        res |= temp;

        res <<= 8;
        temp = bytes[3] & 0xff;
        res |= temp;

        return res;
    }

    inline int64_t bytesToInt64_bigEndian(const unsigned char *b) {
        int64_t temp = 0;
        int64_t res = 0;

        res <<= 8;
        temp = b[0] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[1] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[2] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[3] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[4] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[5] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[6] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[7] & 0xff;
        res |= temp;

        return res;
    }


    inline void int16ToBytes_bigEndian(unsigned char *b, int16_t num) {
        b[0] = (unsigned char) (num >> 8);
        b[1] = (unsigned char) (num);
    }

    inline void int32ToBytes_bigEndian(unsigned char *b, int32_t num) {
        b[0] = (unsigned char) (num >> 24);
        b[1] = (unsigned char) (num >> 16);
        b[2] = (unsigned char) (num >> 8);
        b[3] = (unsigned char) (num);
    }


    inline void int64ToBytes_bigEndian(unsigned char *b, int64_t num) {
        b[0] = (unsigned char) (num >> 56);
        b[1] = (unsigned char) (num >> 48);
        b[2] = (unsigned char) (num >> 40);
        b[3] = (unsigned char) (num >> 32);
        b[4] = (unsigned char) (num >> 24);
        b[5] = (unsigned char) (num >> 16);
        b[6] = (unsigned char) (num >> 8);
        b[7] = (unsigned char) (num);
    }

    std::string floatToBinary(float f) {
        lfloat u;
        u.value = f;
        std::string str(32, '0');
        for (int i = 0; i < 32; i++) {
            str[31 - i] = (u.ivalue % 2) ? '1' : '0';
            u.ivalue >>= 1;
        }
        return str;
    }

    template<class T>
    void truncateArray(T data, size_t n, int byteLen, uchar *&binary) {
        lfloat bytes;
        int b;
        for (size_t i = 0; i < n; i++) {
            bytes.value = data[i];
            for (b = 4 - byteLen; b < 4; b++) {
                *binary++ = bytes.byte[b];
            }
        }
    }

    template<class T>
    void truncateArrayRecover(uchar *binary, size_t n, int byteLen, T *data) {
        lfloat bytes;
        bytes.ivalue = 0;
        int b;
        for (size_t i = 0; i < n; i++) {
            for (b = 4 - byteLen; b < 4; b++) {
                bytes.byte[b] = *binary++;
            }
            data[i] = bytes.value;
        }
    }

    std::vector<uchar> LeadingBitsEncode(float pre, float data) {
        lfloat lfBuf_pre;
        lfloat lfBuf_cur;

        lfBuf_pre.value = pre;
        lfBuf_cur.value = data;
        lfBuf_pre.ivalue = lfBuf_cur.ivalue ^ lfBuf_pre.ivalue;

        std::vector<uchar> bytes;
        int n = 0;
        if (lfBuf_pre.ivalue == 0) {
            n = 0;
        } else if (lfBuf_pre.ivalue >> 8 == 0) {
            n = 1;
        } else if (lfBuf_pre.ivalue >> 16 == 0) {
            n = 2;
        } else if (lfBuf_pre.ivalue >> 24 == 0) {
            n = 3;
        } else {
            n = 4;
        }

        for (int i = 0; i < n; i++) {
            bytes.push_back(lfBuf_cur.byte[i]);
        }
        return bytes;
    }

    float LeadingBitsDecode(float pre, std::vector<uchar> bytes) {
        lfloat lfBuf_pre;
        lfloat lfBuf_cur;

        lfBuf_pre.value = pre;
        lfBuf_cur = lfBuf_pre;

        for (int i = 0; i < bytes.size(); i++) {
            lfBuf_cur.byte[i] = bytes[i];
        }
        return lfBuf_cur.value;
    }

    // *** modified from Sheng's code start ***
    void
    convertIntArray2ByteArray_fast_1b_to_result_sz(const uchar *intArray, size_t intArrayLength, uchar *&compressed_pos)
    {
        size_t byteLength = 0;
        size_t i, j;
        if (intArrayLength % 8 == 0)
            byteLength = intArrayLength / 8;
        else
            byteLength = intArrayLength / 8 + 1;

        size_t n = 0;
        int tmp, type;
        for (i = 0; i < byteLength; i++)
        {
            tmp = 0;
            for (j = 0; j < 8 && n < intArrayLength; j++)
            {
                type = intArray[n];
                if (type == 1)
                    tmp = (tmp | (1 << (7 - j)));
                n++;
            }
            *(compressed_pos++) = (uchar)tmp;
        }
    }

    void convertByteArray2IntArray_fast_1b_sz(size_t intArrayLength, const uchar *&compressed_pos, size_t byteArrayLength, uchar *intArray)
    {
        if (intArrayLength > byteArrayLength * 8)
        {
            printf("Error: intArrayLength > byteArrayLength*8\n");
            printf("intArrayLength=%zu, byteArrayLength = %zu", intArrayLength, byteArrayLength);
            exit(0);
        }
        size_t n = 0, i;
        int tmp;
        for (i = 0; i < byteArrayLength - 1; i++)
        {
            tmp = *(compressed_pos++);
            intArray[n++] = (tmp & 0x80) >> 7;
            intArray[n++] = (tmp & 0x40) >> 6;
            intArray[n++] = (tmp & 0x20) >> 5;
            intArray[n++] = (tmp & 0x10) >> 4;
            intArray[n++] = (tmp & 0x08) >> 3;
            intArray[n++] = (tmp & 0x04) >> 2;
            intArray[n++] = (tmp & 0x02) >> 1;
            intArray[n++] = (tmp & 0x01) >> 0;
        }
        tmp = *(compressed_pos++);
        for (int i = 0; n < intArrayLength; n++, i++)
        {
            intArray[n] = (tmp & (1 << (7 - i))) >> (7 - i);
        }
    }
    // *** modified from Sheng's code end ***


};
#endif //SZ3_BYTEUTIL_HPP
