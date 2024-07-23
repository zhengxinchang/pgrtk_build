impl :: bincode :: Encode for ShmmrSpec
{
    fn encode < __E : :: bincode :: enc :: Encoder >
    (& self, encoder : & mut __E) ->core :: result :: Result < (), :: bincode
    :: error :: EncodeError >
    {
        :: bincode :: Encode :: encode(&self.w, encoder) ?; :: bincode ::
        Encode :: encode(&self.k, encoder) ?; :: bincode :: Encode ::
        encode(&self.r, encoder) ?; :: bincode :: Encode ::
        encode(&self.min_span, encoder) ?; :: bincode :: Encode ::
        encode(&self.sketch, encoder) ?; Ok(())
    }
}