impl :: bincode :: Decode for ShmmrSpec
{
    fn decode < __D : :: bincode :: de :: Decoder > (decoder : & mut __D)
    ->core :: result :: Result < Self, :: bincode :: error :: DecodeError >
    {
        Ok(Self
        {
            w : :: bincode :: Decode :: decode(decoder) ?, k : :: bincode ::
            Decode :: decode(decoder) ?, r : :: bincode :: Decode ::
            decode(decoder) ?, min_span : :: bincode :: Decode ::
            decode(decoder) ?, sketch : :: bincode :: Decode ::
            decode(decoder) ?,
        })
    }
} impl < '__de > :: bincode :: BorrowDecode < '__de > for ShmmrSpec
{
    fn borrow_decode < __D : :: bincode :: de :: BorrowDecoder < '__de > >
    (decoder : & mut __D) ->core :: result :: Result < Self, :: bincode ::
    error :: DecodeError >
    {
        Ok(Self
        {
            w : :: bincode :: BorrowDecode :: borrow_decode(decoder) ?, k : ::
            bincode :: BorrowDecode :: borrow_decode(decoder) ?, r : ::
            bincode :: BorrowDecode :: borrow_decode(decoder) ?, min_span : ::
            bincode :: BorrowDecode :: borrow_decode(decoder) ?, sketch : ::
            bincode :: BorrowDecode :: borrow_decode(decoder) ?,
        })
    }
}