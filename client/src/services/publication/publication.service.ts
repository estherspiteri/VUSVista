import customFetch from "../api/api.service";
import {
  IGetPublicationsByRsidAndWithOptionalTextRequest,
  IGetPublicationsByRsidAndWithOptionalTextResponse,
  IGetPublicationsByVariantIdRequest,
  IGetPublicationsByVariantIdResponse,
} from "./publication.dto";

export class PublicationService {
  async getPublicationsByVariantId(
    input: IGetPublicationsByVariantIdRequest
  ): Promise<IGetPublicationsByVariantIdResponse> {
    const result: IGetPublicationsByVariantIdResponse = await customFetch(
      `/publication/getByVariantId/${input.variantId}`,
      {
        method: "GET",
        headers: {
          "Content-Type": "application/json;charset=UTF-8",
        },
      }
    )
      .then((response) => {
        return response
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async getPublicationsByRsidAndWithOptionalText(
    input: IGetPublicationsByRsidAndWithOptionalTextRequest
  ): Promise<IGetPublicationsByRsidAndWithOptionalTextResponse> {
    const result: IGetPublicationsByVariantIdResponse = await customFetch(
      `/publication/getWithOptionalText/${input.variantId}/${input.rsid}/${input.optionalText}`,
      {
        method: "GET",
        headers: {
          "Content-Type": "application/x-www-form-urlencoded",
        },
      }
    )
      .then((response) => {
        return response
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }
}

export const publicationService = new PublicationService();
