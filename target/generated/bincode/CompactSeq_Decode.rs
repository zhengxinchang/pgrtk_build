impl :: bincode :: Decode for CompactSeq
{
    fn decode < __D : :: bincode :: de :: Decoder > (decoder : & mut __D)
    ->core :: result :: Result < Self, :: bincode :: error :: DecodeError >
    {
        Ok(Self
        {
            source : :: bincode :: Decode :: decode(decoder) ?, name : ::
            bincode :: Decode :: decode(decoder) ?, id : :: bincode :: Decode
            :: decode(decoder) ?, seq_frag_range : :: bincode :: Decode ::
            decode(decoder) ?, len : :: bincode :: Decode :: decode(decoder)
            ?,
        })
    }
} impl < '__de > :: bincode :: BorrowDecode < '__de > for CompactSeq
{
    fn borrow_decode < __D : :: bincode :: de :: BorrowDecoder < '__de > >
    (decoder : & mut __D) ->core :: result :: Result < Self, :: bincode ::
    error :: DecodeError >
    {
        Ok(Self
        {
            source : :: bincode :: BorrowDecode :: borrow_decode(decoder) ?,
            name : :: bincode :: BorrowDecode :: borrow_decode(decoder) ?, id
            : :: bincode :: BorrowDecode :: borrow_decode(decoder) ?,
            seq_frag_range : :: bincode :: BorrowDecode ::
            borrow_decode(decoder) ?, len : :: bincode :: BorrowDecode ::
            borrow_decode(decoder) ?,
        })
    }
}