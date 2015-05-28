# steganography
Hiding of encrypted text in images using DCT
1) Text is encrypted usig standard RCA algorithm 
2) DCT of image is taken and quantisation is done 
3) The text is then embedded
4) Huffman coding is then used for decreasing the length of the bits to be stored
5) This is similar to jpeg compression

Decoding is done by just reversing the above process
