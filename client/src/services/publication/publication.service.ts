import {
  IGetPublicationsByVariantIdRequest,
  IGetPublicationsByVariantIdResponse,
} from "./publication.dto";

export class PublicationService {
  async getPublicationsByVariantId(
    input: IGetPublicationsByVariantIdRequest
  ): Promise<IGetPublicationsByVariantIdResponse> {
    const result: IGetPublicationsByVariantIdResponse = await fetch(
      `/publication/getByVariantId/${input.variantId}`,
      {
        method: "GET",
        headers: {
          "Content-Type": "application/json;charset=UTF-8",
        },
      }
    )
      .then((response: Response) => {
        return response.json();
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }
}

export const publicationService = new PublicationService();
